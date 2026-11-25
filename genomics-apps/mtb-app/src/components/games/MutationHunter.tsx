'use client';

import { useCallback, useEffect, useMemo, useState } from 'react';
import { CODON_TABLE, DNA_BASES, SENSE_CODONS, THREE_LETTER, translate } from './geneticCode';
import { readBest, saveBest } from './bestScore';

// Mutation Hunter: compare a patient's coding sequence with the reference and name the kind of variant

type MutationClass = 'Silent' | 'Missense' | 'Nonsense' | 'Frameshift' | 'In-frame deletion';
const CLASSES: MutationClass[] = ['Silent', 'Missense', 'Nonsense', 'Frameshift', 'In-frame deletion'];

const CODON_COUNT = 8;
const ROUNDS = 10;
const ROUND_SECONDS = 30;

// One column of the alignment; a gap is '-'
interface AlignedBase { reference: string; patient: string; changed: boolean }

interface Puzzle {
  kind: MutationClass;
  alignment: AlignedBase[];
  referenceProtein: string[];
  patientProtein: string[];
  hgvs: string;
}

const pick = <T,>(items: readonly T[]): T => items[Math.floor(Math.random() * items.length)];

function randomReference(): string {
  return 'ATG' + Array.from({ length: CODON_COUNT - 1 }, () => pick(SENSE_CODONS)).join('');
}

function substitutionPuzzle(kind: MutationClass): Puzzle | null {
  for (let attempt = 0; attempt < 200; attempt++) {
    const reference = randomReference();
    const position = 3 + Math.floor(Math.random() * (reference.length - 3));
    const base = pick(DNA_BASES.filter((b) => b !== reference[position]));
    const patient = reference.slice(0, position) + base + reference.slice(position + 1);
    const codonIndex = Math.floor(position / 3);
    const before = CODON_TABLE[reference.slice(codonIndex * 3, codonIndex * 3 + 3)];
    const after = CODON_TABLE[patient.slice(codonIndex * 3, codonIndex * 3 + 3)];
    const actual: MutationClass = after === before ? 'Silent' : after === '*' ? 'Nonsense' : 'Missense';
    if (actual !== kind) continue;
    const residue = `${THREE_LETTER[before]}${codonIndex + 1}`;
    return {
      kind,
      alignment: [...reference].map((b, i) => ({ reference: b, patient: patient[i], changed: i === position })),
      referenceProtein: translate(reference),
      patientProtein: translate(patient),
      hgvs: `p.${residue}${kind === 'Silent' ? '=' : THREE_LETTER[after]}`,
    };
  }
  return null;
}

function frameshiftPuzzle(): Puzzle {
  const reference = randomReference();
  const position = 3 + Math.floor(Math.random() * (reference.length - 6));
  const insertion = Math.random() < 0.5;
  const alignment: AlignedBase[] = [];
  let patient = '';
  if (insertion) {
    const base = pick(DNA_BASES);
    patient = reference.slice(0, position) + base + reference.slice(position);
    [...reference].forEach((b, i) => {
      if (i === position) alignment.push({ reference: '-', patient: base, changed: true });
      alignment.push({ reference: b, patient: b, changed: false });
    });
  } else {
    patient = reference.slice(0, position) + reference.slice(position + 1);
    [...reference].forEach((b, i) => alignment.push({ reference: b, patient: i === position ? '-' : b, changed: i === position }));
  }
  const codonIndex = Math.floor(position / 3);
  return {
    kind: 'Frameshift',
    alignment,
    referenceProtein: translate(reference),
    patientProtein: translate(patient),
    hgvs: `p.${THREE_LETTER[CODON_TABLE[reference.slice(codonIndex * 3, codonIndex * 3 + 3)]]}${codonIndex + 1}fs`,
  };
}

function inFrameDeletionPuzzle(): Puzzle {
  const reference = randomReference();
  const codonIndex = 1 + Math.floor(Math.random() * (CODON_COUNT - 2));
  const start = codonIndex * 3;
  const patient = reference.slice(0, start) + reference.slice(start + 3);
  return {
    kind: 'In-frame deletion',
    alignment: [...reference].map((b, i) => {
      const deleted = i >= start && i < start + 3;
      return { reference: b, patient: deleted ? '-' : b, changed: deleted };
    }),
    referenceProtein: translate(reference),
    patientProtein: translate(patient),
    hgvs: `p.${THREE_LETTER[CODON_TABLE[reference.slice(start, start + 3)]]}${codonIndex + 1}del`,
  };
}

function makePuzzle(): Puzzle {
  const kind = pick(CLASSES);
  if (kind === 'Frameshift') return frameshiftPuzzle();
  if (kind === 'In-frame deletion') return inFrameDeletionPuzzle();
  return substitutionPuzzle(kind) ?? frameshiftPuzzle();
}

const EXPLANATIONS: Record<MutationClass, string> = {
  Silent: 'The new codon codes for the same amino acid, so the protein is unchanged.',
  Missense: 'The new codon codes for a different amino acid: one residue of the protein changes.',
  Nonsense: 'The new codon is a stop codon, so the protein is cut short.',
  Frameshift: 'One base gained or lost shifts the reading frame; every codon after it is read wrongly, usually until an early stop.',
  'In-frame deletion': 'Three bases (one whole codon) are lost: the protein loses one amino acid but keeps its frame.',
};

function SequenceRow({ label, bases, revealed }: { label: string; bases: { base: string; changed: boolean }[]; revealed: boolean }) {
  return (
    <div className="flex items-center gap-3">
      <span className="w-20 shrink-0 text-right text-xs font-semibold uppercase tracking-wide text-slate-500">{label}</span>
      <div className="flex flex-wrap font-mono text-lg">
        {bases.map(({ base, changed }, i) => (
          <span
            key={i}
            className={`inline-block w-6 text-center ${i % 3 === 2 ? 'mr-2' : ''} ${changed ? (revealed ? 'rounded bg-rose-500 font-bold text-white' : 'rounded bg-amber-200 font-bold') : 'text-slate-800'}`}
          >
            {base}
          </span>
        ))}
      </div>
    </div>
  );
}

function ProteinRow({ label, residues, compareTo }: { label: string; residues: string[]; compareTo?: string[] }) {
  return (
    <div className="flex items-center gap-3">
      <span className="w-20 shrink-0 text-right text-xs font-semibold uppercase tracking-wide text-slate-500">{label}</span>
      <div className="flex flex-wrap gap-1 font-mono text-sm">
        {residues.map((residue, i) => {
          const differs = compareTo && compareTo[i] !== residue;
          return (
            <span key={i} className={`rounded px-1.5 py-0.5 ${residue === '*' ? 'bg-slate-900 text-white' : differs ? 'bg-rose-100 text-rose-800' : 'bg-slate-100 text-slate-700'}`}>
              {THREE_LETTER[residue]}
            </span>
          );
        })}
      </div>
    </div>
  );
}

function CodonTable() {
  const groups = useMemo(() => {
    const byAmino: Record<string, string[]> = {};
    Object.entries(CODON_TABLE).forEach(([codon, amino]) => (byAmino[amino] ??= []).push(codon));
    return Object.entries(byAmino).sort(([a], [b]) => THREE_LETTER[a].localeCompare(THREE_LETTER[b]));
  }, []);
  return (
    <div className="grid grid-cols-2 gap-x-4 gap-y-1 text-xs sm:grid-cols-3 lg:grid-cols-4">
      {groups.map(([amino, codons]) => (
        <div key={amino} className="flex gap-2">
          <span className={`w-9 font-semibold ${amino === '*' ? 'text-rose-700' : 'text-slate-800'}`}>{THREE_LETTER[amino]}</span>
          <span className="font-mono text-slate-600">{codons.join(' ')}</span>
        </div>
      ))}
    </div>
  );
}

export default function MutationHunter() {
  const [round, setRound] = useState(0);
  const [puzzle, setPuzzle] = useState<Puzzle | null>(null);
  const [answer, setAnswer] = useState<MutationClass | 'timeout' | null>(null);
  const [score, setScore] = useState(0);
  const [streak, setStreak] = useState(0);
  const [secondsLeft, setSecondsLeft] = useState(ROUND_SECONDS);
  const [showTable, setShowTable] = useState(false);
  const [best, setBest] = useState(0);
  const [finished, setFinished] = useState(false);
  const [newBest, setNewBest] = useState(false);

  useEffect(() => setBest(readBest('mutation-hunter')), []);

  const start = useCallback(() => {
    setRound(1);
    setScore(0);
    setStreak(0);
    setFinished(false);
    setAnswer(null);
    setSecondsLeft(ROUND_SECONDS);
    setPuzzle(makePuzzle());
  }, []);

  const submit = useCallback((choice: MutationClass | 'timeout') => {
    if (!puzzle || answer) return;
    const correct = choice === puzzle.kind;
    setAnswer(choice);
    if (correct) {
      const nextStreak = streak + 1;
      setStreak(nextStreak);
      setScore((s) => s + 100 + secondsLeft * 5 + (nextStreak - 1) * 25);
    } else {
      setStreak(0);
    }
  }, [puzzle, answer, streak, secondsLeft]);

  // Count down while the round is open
  useEffect(() => {
    if (!puzzle || answer || finished) return;
    if (secondsLeft <= 0) {
      submit('timeout');
      return;
    }
    const timer = setTimeout(() => setSecondsLeft((s) => s - 1), 1000);
    return () => clearTimeout(timer);
  }, [puzzle, answer, finished, secondsLeft, submit]);

  const next = () => {
    if (round >= ROUNDS) {
      setFinished(true);
      setNewBest(score > best);
      if (score > best) {
        setBest(score);
        saveBest('mutation-hunter', score);
      }
      return;
    }
    setRound((r) => r + 1);
    setAnswer(null);
    setSecondsLeft(ROUND_SECONDS);
    setPuzzle(makePuzzle());
  };

  if (!puzzle || finished) {
    return (
      <div className="space-y-4 text-center">
        {finished ? (
          <>
            <p className="text-2xl font-bold text-slate-900">Final score: {score}</p>
            <p className="text-sm text-slate-600">{newBest ? 'A new personal best!' : `Your best: ${best}`}</p>
          </>
        ) : (
          <p className="mx-auto max-w-xl text-sm leading-6 text-slate-600">
            Each round shows a short stretch of a gene&apos;s coding sequence, the reference above and the patient below,
            with the variant highlighted. Name the kind of variant before the clock runs out. Faster answers and streaks
            score more. {ROUNDS} rounds. {best > 0 && <>Your best: <strong>{best}</strong>.</>}
          </p>
        )}
        <button type="button" onClick={start} className="rounded-lg bg-slate-900 px-6 py-2.5 text-sm font-semibold text-white hover:bg-slate-700">
          {finished ? 'Play again' : 'Start'}
        </button>
      </div>
    );
  }

  const revealed = answer !== null;
  const correct = answer === puzzle.kind;
  return (
    <div className="space-y-5">
      <div className="flex flex-wrap items-center justify-between gap-3 text-sm">
        <span className="font-semibold text-slate-700">Round {round} / {ROUNDS}</span>
        <span className="text-slate-600">Score <strong className="text-slate-900">{score}</strong>{streak > 1 && <span className="ml-2 rounded-full bg-amber-100 px-2 py-0.5 text-xs font-semibold text-amber-800">🔥 {streak} streak</span>}</span>
        <div className="h-2 w-40 overflow-hidden rounded-full bg-slate-200" aria-label={`${secondsLeft} seconds left`}>
          <div className={`h-full transition-all ${secondsLeft <= 8 ? 'bg-rose-500' : 'bg-emerald-500'}`} style={{ width: `${(secondsLeft / ROUND_SECONDS) * 100}%` }} />
        </div>
      </div>

      <div className="space-y-2 overflow-x-auto rounded-xl bg-white p-4 ring-1 ring-slate-200">
        <SequenceRow label="Reference" bases={puzzle.alignment.map((a) => ({ base: a.reference, changed: a.changed }))} revealed={revealed} />
        <SequenceRow label="Patient" bases={puzzle.alignment.map((a) => ({ base: a.patient, changed: a.changed }))} revealed={revealed} />
        {revealed && (
          <div className="mt-3 space-y-2 border-t border-slate-100 pt-3">
            <ProteinRow label="Ref protein" residues={puzzle.referenceProtein} />
            <ProteinRow label="Patient" residues={puzzle.patientProtein} compareTo={puzzle.referenceProtein} />
          </div>
        )}
      </div>

      <div className="flex flex-wrap gap-2">
        {CLASSES.map((kind) => {
          const isAnswer = revealed && kind === puzzle.kind;
          const isWrongPick = revealed && kind === answer && !correct;
          return (
            <button
              key={kind}
              type="button"
              disabled={revealed}
              onClick={() => submit(kind)}
              className={`rounded-lg border px-4 py-2 text-sm font-semibold transition-colors ${isAnswer ? 'border-emerald-600 bg-emerald-600 text-white' : isWrongPick ? 'border-rose-500 bg-rose-500 text-white' : 'border-slate-300 bg-white text-slate-700 enabled:hover:bg-slate-50 disabled:opacity-60'}`}
            >
              {kind}
            </button>
          );
        })}
      </div>

      {revealed && (
        <div className={`flex flex-wrap items-start justify-between gap-3 rounded-xl p-4 text-sm ${correct ? 'bg-emerald-50 text-emerald-900' : 'bg-rose-50 text-rose-900'}`}>
          <div className="space-y-1">
            <p className="font-semibold">
              {correct ? 'Correct!' : answer !== 'timeout' ? `Not quite — it's ${puzzle.kind.toLowerCase()}.` : `Time's up — it's ${puzzle.kind.toLowerCase()}.`}{' '}
              <span className="font-mono font-normal">{puzzle.hgvs}</span>
            </p>
            <p>{EXPLANATIONS[puzzle.kind]}</p>
          </div>
          <button type="button" onClick={next} className="rounded-lg bg-slate-900 px-4 py-2 font-semibold text-white hover:bg-slate-700">
            {round >= ROUNDS ? 'See score' : 'Next'}
          </button>
        </div>
      )}

      <div>
        <button type="button" onClick={() => setShowTable((s) => !s)} className="text-sm font-medium text-blue-600 hover:underline">
          {showTable ? 'Hide' : 'Show'} the genetic code
        </button>
        {showTable && <div className="mt-3 rounded-xl bg-white p-4 ring-1 ring-slate-200"><CodonTable /></div>}
      </div>
    </div>
  );
}
