'use client';

import { useState } from 'react';
import { Dna, LucideIcon, Microscope, Pill } from 'lucide-react';
import MutationHunter from './MutationHunter';
import TargetMatch from './TargetMatch';
import CloneControl from './CloneControl';

type GameId = 'mutation-hunter' | 'target-match' | 'clone-control';

const GAMES: { id: GameId; title: string; tagline: string; icon: LucideIcon; accent: string }[] = [
  { id: 'mutation-hunter', title: 'Mutation Hunter', tagline: 'Spot the variant and name its effect on the protein', icon: Dna, accent: 'from-sky-500 to-cyan-400' },
  { id: 'target-match', title: 'Target Match', tagline: 'A memory game pairing alterations with targeted therapies', icon: Pill, accent: 'from-violet-500 to-fuchsia-400' },
  { id: 'clone-control', title: 'Clone Control', tagline: 'Outwit tumor evolution: dose drugs before resistance takes over', icon: Microscope, accent: 'from-rose-500 to-amber-400' },
];

// Games for a break between cases; they don't use the selected patient
export default function GamesView() {
  const [game, setGame] = useState<GameId>('mutation-hunter');

  return (
    <div className="space-y-5">
      <div className="grid gap-3 sm:grid-cols-3">
        {GAMES.map(({ id, title, tagline, icon: Icon, accent }) => {
          const isActive = game === id;
          return (
            <button
              key={id}
              type="button"
              onClick={() => setGame(id)}
              aria-pressed={isActive}
              className={`flex items-start gap-3 rounded-2xl bg-white p-4 text-left shadow-sm transition-all ${isActive ? 'ring-2 ring-slate-900' : 'ring-1 ring-slate-200 hover:ring-slate-400'}`}
            >
              <span className={`flex h-10 w-10 shrink-0 items-center justify-center rounded-xl bg-gradient-to-br text-white ${accent}`}>
                <Icon size={20} aria-hidden="true" />
              </span>
              <span>
                <span className="block font-semibold text-slate-900">{title}</span>
                <span className="block text-xs leading-5 text-slate-500">{tagline}</span>
              </span>
            </button>
          );
        })}
      </div>

      <section className="rounded-2xl bg-white p-5 shadow-sm ring-1 ring-slate-200">
        <h2 className="mb-4 text-lg font-bold text-slate-900">{GAMES.find((g) => g.id === game)?.title}</h2>
        {game === 'mutation-hunter' && <MutationHunter />}
        {game === 'target-match' && <TargetMatch />}
        {game === 'clone-control' && <CloneControl />}
      </section>

      <p className="text-center text-xs text-slate-400">For fun and learning only — not clinical guidance.</p>
    </div>
  );
}
