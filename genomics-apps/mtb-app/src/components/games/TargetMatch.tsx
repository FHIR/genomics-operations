'use client';

import { useEffect, useState } from 'react';
import { readBest, saveBest } from './bestScore';

// Target Match: a memory game pairing each genomic alteration with a targeted therapy approved for it.
// Each drug appears once, so every card has exactly one partner.

interface TargetPair { alteration: string; drug: string; note: string }

const PAIRS: TargetPair[] = [
  { alteration: 'EGFR L858R', drug: 'Osimertinib', note: 'Third-generation EGFR inhibitor, first-line in EGFR-mutant NSCLC.' },
  { alteration: 'EGFR exon 20 insertion', drug: 'Amivantamab', note: 'EGFR–MET bispecific antibody; exon 20 insertions resist most EGFR TKIs.' },
  { alteration: 'BRAF V600E', drug: 'Dabrafenib + trametinib', note: 'BRAF plus MEK inhibition, used in melanoma, NSCLC and other V600E tumors.' },
  { alteration: 'KRAS G12C', drug: 'Sotorasib', note: 'Locks KRAS G12C in its inactive state; the first approved KRAS inhibitor.' },
  { alteration: 'ALK fusion', drug: 'Alectinib', note: 'CNS-penetrant ALK inhibitor for ALK-rearranged NSCLC.' },
  { alteration: 'ROS1 fusion', drug: 'Crizotinib', note: 'Multi-kinase inhibitor active against ROS1-rearranged NSCLC.' },
  { alteration: 'NTRK fusion', drug: 'Larotrectinib', note: 'Tumor-agnostic approval: works on NTRK fusions in any cancer type.' },
  { alteration: 'RET fusion', drug: 'Selpercatinib', note: 'Selective RET inhibitor for RET-altered lung and thyroid cancers.' },
  { alteration: 'MET exon 14 skipping', drug: 'Capmatinib', note: 'MET inhibitor for NSCLC with MET exon 14 skipping.' },
  { alteration: 'ERBB2 (HER2) amplification', drug: 'Trastuzumab', note: 'Anti-HER2 antibody; transformed HER2-positive breast cancer care.' },
  { alteration: 'BCR::ABL1 fusion', drug: 'Imatinib', note: 'The archetypal targeted therapy, for Philadelphia chromosome-positive CML.' },
  { alteration: 'BRCA1/2 loss', drug: 'Olaparib', note: 'PARP inhibitor exploiting synthetic lethality in homologous-recombination–deficient tumors.' },
  { alteration: 'IDH1 R132', drug: 'Ivosidenib', note: 'IDH1 inhibitor for IDH1-mutant AML and cholangiocarcinoma.' },
  { alteration: 'IDH2 R140Q', drug: 'Enasidenib', note: 'IDH2 inhibitor that lets IDH2-mutant AML blasts differentiate.' },
  { alteration: 'FLT3-ITD', drug: 'Midostaurin', note: 'Added to induction chemotherapy for FLT3-mutant AML.' },
  { alteration: 'PIK3CA H1047R', drug: 'Alpelisib', note: 'PI3Kα inhibitor for PIK3CA-mutant HR-positive breast cancer.' },
  { alteration: 'FGFR2 fusion', drug: 'Pemigatinib', note: 'FGFR inhibitor for FGFR2-rearranged cholangiocarcinoma.' },
  { alteration: 'MSI-high / dMMR', drug: 'Pembrolizumab', note: 'Checkpoint inhibitor; hypermutated MSI-high tumors are highly immunogenic.' },
];

const PAIRS_PER_GAME = 8;

interface Card { id: number; pairIndex: number; side: 'alteration' | 'drug'; text: string }

function newDeck(): { cards: Card[]; pairs: TargetPair[] } {
  const pairs = [...PAIRS].sort(() => Math.random() - 0.5).slice(0, PAIRS_PER_GAME);
  const cards = pairs.flatMap((pair, pairIndex) => [
    { pairIndex, side: 'alteration' as const, text: pair.alteration },
    { pairIndex, side: 'drug' as const, text: pair.drug },
  ]);
  for (let i = cards.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [cards[i], cards[j]] = [cards[j], cards[i]];
  }
  return { cards: cards.map((card, id) => ({ ...card, id })), pairs };
}

export default function TargetMatch() {
  const [deck, setDeck] = useState<{ cards: Card[]; pairs: TargetPair[] } | null>(null);
  const [flipped, setFlipped] = useState<number[]>([]);
  const [matched, setMatched] = useState<number[]>([]);
  const [moves, setMoves] = useState(0);
  const [seconds, setSeconds] = useState(0);
  const [best, setBest] = useState(0);

  // Deal on the client only, so the shuffle doesn't differ between server and browser renders
  useEffect(() => {
    setDeck(newDeck());
    setBest(readBest('target-match'));
  }, []);

  const done = deck !== null && matched.length === deck.pairs.length;
  const score = Math.max(0, 1000 - (moves - PAIRS_PER_GAME) * 40 - seconds * 2);

  useEffect(() => {
    if (!deck || done || moves === 0) return;
    const timer = setInterval(() => setSeconds((s) => s + 1), 1000);
    return () => clearInterval(timer);
  }, [deck, done, moves]);

  useEffect(() => {
    if (done && score > best) {
      setBest(score);
      saveBest('target-match', score);
    }
  }, [done, score, best]);

  const reset = () => {
    setDeck(newDeck());
    setFlipped([]);
    setMatched([]);
    setMoves(0);
    setSeconds(0);
  };

  const flip = (card: Card) => {
    if (!deck || flipped.length === 2 || flipped.includes(card.id) || matched.includes(card.pairIndex)) return;
    const next = [...flipped, card.id];
    setFlipped(next);
    if (next.length < 2) return;
    setMoves((m) => m + 1);
    const [first, second] = next.map((id) => deck.cards[id]);
    if (first.pairIndex === second.pairIndex) {
      setMatched((m) => [...m, first.pairIndex]);
      setFlipped([]);
    } else {
      setTimeout(() => setFlipped([]), 1100);
    }
  };

  if (!deck) return null;
  const lastMatch = matched.length ? deck.pairs[matched[matched.length - 1]] : null;

  return (
    <div className="space-y-5">
      <div className="flex flex-wrap items-center justify-between gap-3 text-sm text-slate-600">
        <span>Match each alteration with its targeted therapy.</span>
        <span className="flex gap-4">
          <span>Moves <strong className="text-slate-900">{moves}</strong></span>
          <span>Time <strong className="text-slate-900">{seconds}s</strong></span>
          <span>Pairs <strong className="text-slate-900">{matched.length}/{deck.pairs.length}</strong></span>
        </span>
      </div>

      <div className="grid grid-cols-2 gap-3 sm:grid-cols-4">
        {deck.cards.map((card) => {
          const isMatched = matched.includes(card.pairIndex);
          const isUp = isMatched || flipped.includes(card.id);
          return (
            <button
              key={card.id}
              type="button"
              onClick={() => flip(card)}
              aria-label={isUp ? card.text : 'Hidden card'}
              className={`flex h-24 items-center justify-center rounded-xl p-3 text-center text-sm font-semibold shadow-sm transition-all duration-300 ${
                isMatched
                  ? 'bg-emerald-100 text-emerald-900 ring-2 ring-emerald-400'
                  : isUp
                    ? card.side === 'alteration' ? 'bg-sky-50 text-sky-900 ring-2 ring-sky-400' : 'bg-violet-50 text-violet-900 ring-2 ring-violet-400'
                    : 'bg-gradient-to-br from-slate-700 to-slate-900 text-slate-400 hover:from-slate-600 hover:to-slate-800'
              }`}
            >
              {isUp ? (
                <span>
                  <span className="mb-1 block text-[10px] font-bold uppercase tracking-wider opacity-60">{card.side === 'alteration' ? 'Alteration' : 'Therapy'}</span>
                  {card.side === 'alteration' ? <span className="font-mono">{card.text}</span> : card.text}
                </span>
              ) : (
                <span className="text-2xl" aria-hidden="true">🧬</span>
              )}
            </button>
          );
        })}
      </div>

      {done ? (
        <div className="flex flex-wrap items-center justify-between gap-3 rounded-xl bg-emerald-50 p-4 text-sm text-emerald-900">
          <span><strong>Board cleared!</strong> {moves} moves in {seconds}s — score {score}{score >= best ? ' (personal best)' : ` (best ${best})`}.</span>
          <button type="button" onClick={reset} className="rounded-lg bg-slate-900 px-4 py-2 font-semibold text-white hover:bg-slate-700">Deal again</button>
        </div>
      ) : lastMatch && (
        <p className="rounded-xl bg-white p-4 text-sm text-slate-700 ring-1 ring-slate-200">
          <strong className="font-mono">{lastMatch.alteration}</strong> → <strong>{lastMatch.drug}</strong>: {lastMatch.note}
        </p>
      )}
    </div>
  );
}
