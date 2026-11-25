'use client';

import { useCallback, useEffect, useRef, useState } from 'react';
import { CloneState, countClones, createCloneState, GRID_HEIGHT, GRID_WIDTH, MAX_BURDEN, stepCloneState, Treatment, WIN_WEEKS } from './cloneControlEngine';
import { readBest, saveBest } from './bestScore';

const CELL = 14;
// Empty, sensitive, resistant to A, resistant to B, resistant to both
const CELL_COLORS = ['#f8fafc', '#38bdf8', '#f59e0b', '#a855f7', '#e11d48'];

const TREATMENTS: { value: Treatment; label: string; key: string; hint: string }[] = [
  { value: 'A', label: 'Drug A', key: '1', hint: 'Kills blue and purple cells' },
  { value: 'B', label: 'Drug B', key: '2', hint: 'Kills blue and orange cells' },
  { value: 'AB', label: 'A + B', key: '3', hint: 'Kills all but red; very toxic' },
  { value: 'holiday', label: 'Drug holiday', key: '4', hint: 'Patient recovers; tumor regrows' },
];

const SPEEDS = [{ label: 'Slow', ms: 700 }, { label: 'Normal', ms: 400 }, { label: 'Fast', ms: 180 }];

const OUTCOME_TEXT: Record<Exclude<CloneState['outcome'], 'playing'>, { title: string; body: string; good: boolean }> = {
  survived: { title: `Two years of disease control!`, body: 'You kept the tumor in check without exhausting the patient. Sparing sensitive cells keeps resistant clones crowded out — the idea behind adaptive therapy.', good: true },
  cleared: { title: 'Complete response!', body: 'Every tumor cell is gone. Rare and wonderful.', good: true },
  progression: { title: 'Disease progression', body: 'The tumor outgrew its space. Resistant clones thrive when continuous treatment wipes out their sensitive competitors.', good: false },
  toxicity: { title: 'Treatment-limiting toxicity', body: 'The patient could no longer tolerate treatment. Drug holidays let them recover.', good: false },
};

function Meter({ label, value, max, danger }: { label: string; value: number; max: number; danger: number }) {
  const fraction = Math.min(1, value / max);
  return (
    <div className="min-w-[140px] flex-1">
      <div className="mb-1 flex justify-between text-xs font-semibold text-slate-600"><span>{label}</span><span>{Math.round(fraction * 100)}%</span></div>
      <div className="h-2.5 overflow-hidden rounded-full bg-slate-200">
        <div className={`h-full transition-all ${fraction >= danger ? 'bg-rose-500' : fraction >= danger * 0.7 ? 'bg-amber-500' : 'bg-emerald-500'}`} style={{ width: `${fraction * 100}%` }} />
      </div>
    </div>
  );
}

export default function CloneControl() {
  const [state, setState] = useState<CloneState | null>(null);
  const [treatment, setTreatment] = useState<Treatment>('holiday');
  const [running, setRunning] = useState(false);
  const [speed, setSpeed] = useState(1);
  const [best, setBest] = useState(0);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const treatmentRef = useRef(treatment);
  treatmentRef.current = treatment;

  useEffect(() => {
    setState(createCloneState());
    setBest(readBest('clone-control'));
  }, []);

  const restart = () => {
    setState(createCloneState());
    setTreatment('holiday');
    setRunning(false);
  };

  // Advance one week per tick while running
  useEffect(() => {
    if (!running) return;
    const timer = setInterval(() => {
      setState((current) => (current ? stepCloneState(current, treatmentRef.current) : current));
    }, SPEEDS[speed].ms);
    return () => clearInterval(timer);
  }, [running, speed]);

  useEffect(() => {
    if (state && state.outcome !== 'playing') {
      setRunning(false);
      const weeks = state.outcome === 'cleared' ? WIN_WEEKS : state.week;
      if (weeks > best) {
        setBest(weeks);
        saveBest('clone-control', weeks);
      }
    }
  }, [state, best]);

  useEffect(() => {
    const context = canvasRef.current?.getContext('2d');
    if (!context || !state) return;
    context.fillStyle = CELL_COLORS[0];
    context.fillRect(0, 0, GRID_WIDTH * CELL, GRID_HEIGHT * CELL);
    state.cells.forEach((value, index) => {
      if (!value) return;
      context.fillStyle = CELL_COLORS[value];
      context.beginPath();
      context.arc((index % GRID_WIDTH) * CELL + CELL / 2, Math.floor(index / GRID_WIDTH) * CELL + CELL / 2, CELL / 2 - 1, 0, Math.PI * 2);
      context.fill();
    });
  }, [state]);

  // Number keys pick the treatment; space pauses
  const onKey = useCallback((event: KeyboardEvent) => {
    const target = event.target as HTMLElement | null;
    if (target && ['INPUT', 'TEXTAREA', 'SELECT'].includes(target.tagName)) return;
    const choice = TREATMENTS.find((t) => t.key === event.key);
    if (choice) setTreatment(choice.value);
    if (event.key === ' ' && state?.outcome === 'playing') {
      event.preventDefault();
      setRunning((r) => !r);
    }
  }, [state]);
  useEffect(() => {
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [onKey]);

  if (!state) return null;
  const counts = countClones(state.cells);
  const outcome = state.outcome === 'playing' ? null : OUTCOME_TEXT[state.outcome];

  return (
    <div className="space-y-4">
      <p className="text-sm leading-6 text-slate-600">
        Keep the patient alive for {WIN_WEEKS} weeks. Each week, choose a treatment. Drugs kill the cells that are sensitive to them,
        but every division risks a resistance mutation — and resistant cells grow a little more slowly than sensitive ones.
        Don&apos;t let the tumor fill {Math.round(MAX_BURDEN * 100)}% of the field, or toxicity reach 100%. Keys 1–4 choose a treatment; space pauses.
      </p>

      <div className="flex flex-wrap items-end gap-4">
        <div className="text-sm"><span className="text-slate-500">Week</span> <strong className="text-2xl text-slate-900">{state.week}</strong><span className="text-slate-400"> / {WIN_WEEKS}</span></div>
        <Meter label="Tumor burden" value={counts.total / state.cells.length} max={MAX_BURDEN} danger={0.85} />
        <Meter label="Toxicity" value={state.toxicity} max={100} danger={0.85} />
        {best > 0 && <span className="text-xs text-slate-500">Best: {best >= WIN_WEEKS ? 'won' : `${best} weeks`}</span>}
      </div>

      <div className="relative">
        <canvas
          ref={canvasRef}
          width={GRID_WIDTH * CELL}
          height={GRID_HEIGHT * CELL}
          className="h-auto w-full max-w-[672px] rounded-xl bg-slate-50 ring-1 ring-slate-200"
        />
        {outcome && (
          <div className="absolute inset-0 flex max-w-[672px] items-center justify-center rounded-xl bg-white/80 p-6 backdrop-blur-sm">
            <div className="max-w-md space-y-3 text-center">
              <p className={`text-xl font-bold ${outcome.good ? 'text-emerald-700' : 'text-rose-700'}`}>{outcome.title}</p>
              <p className="text-sm text-slate-700">{outcome.body}</p>
              <p className="text-sm text-slate-500">Weeks survived: {state.week}</p>
              <button type="button" onClick={restart} className="rounded-lg bg-slate-900 px-5 py-2 text-sm font-semibold text-white hover:bg-slate-700">New patient</button>
            </div>
          </div>
        )}
      </div>

      <div className="flex flex-wrap gap-x-4 gap-y-1 text-xs text-slate-600">
        {[['Sensitive', counts.sensitive, 1], ['Resistant to A', counts.resistantA, 2], ['Resistant to B', counts.resistantB, 3], ['Resistant to both', counts.resistantBoth, 4]].map(([label, count, color]) => (
          <span key={label as string} className="inline-flex items-center gap-1.5">
            <span className="h-3 w-3 rounded-full" style={{ background: CELL_COLORS[color as number] }} />
            {label} <strong className="text-slate-800">{count}</strong>
          </span>
        ))}
      </div>

      <div className="flex flex-wrap items-center gap-2">
        {TREATMENTS.map(({ value, label, key, hint }) => (
          <button
            key={value}
            type="button"
            title={hint}
            onClick={() => setTreatment(value)}
            className={`rounded-lg border px-4 py-2 text-sm font-semibold ${treatment === value ? 'border-slate-900 bg-slate-900 text-white' : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-50'}`}
          >
            <span className="mr-1.5 text-xs opacity-60">{key}</span>{label}
          </button>
        ))}
        <span className="mx-2 h-6 w-px bg-slate-300" aria-hidden="true" />
        <button
          type="button"
          disabled={state.outcome !== 'playing'}
          onClick={() => setRunning((r) => !r)}
          className="rounded-lg bg-emerald-600 px-5 py-2 text-sm font-semibold text-white hover:bg-emerald-500 disabled:opacity-50"
        >
          {running ? 'Pause' : state.week === 0 ? 'Start' : 'Resume'}
        </button>
        <select
          value={speed}
          onChange={(event) => setSpeed(Number(event.target.value))}
          aria-label="Speed"
          className="rounded-lg border border-slate-300 bg-white px-2 py-2 text-sm"
        >
          {SPEEDS.map((s, i) => <option key={s.label} value={i}>{s.label}</option>)}
        </select>
      </div>
    </div>
  );
}
