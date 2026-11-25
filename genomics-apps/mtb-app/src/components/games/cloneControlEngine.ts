// Clone Control: a tumor grows on a grid of cells; the player chooses each week's treatment.
// Each cell is 0 (empty) or 1 + resistance bits (bit 1 = resistant to drug A, bit 2 = resistant to drug B).
// Resistance carries a fitness cost, so sensitive cells out-compete resistant ones for space when the drugs stop,
// which is the idea behind adaptive therapy.

export type Treatment = 'A' | 'B' | 'AB' | 'holiday';

export const GRID_WIDTH = 48;
export const GRID_HEIGHT = 30;
export const WIN_WEEKS = 104;
export const MAX_BURDEN = 0.5;

const DIVISION_RATE = [0, 0.16, 0.1, 0.1, 0.07];
const MUTATION_RATE = 0.0015;
const KILL_RATE = 0.12;
const TOXICITY_PER_WEEK: Record<Treatment, number> = { A: 1.2, B: 1.2, AB: 4, holiday: -4 };

export interface CloneState {
  cells: Uint8Array;
  week: number;
  toxicity: number;
  outcome: 'playing' | 'progression' | 'toxicity' | 'survived' | 'cleared';
}

export interface CloneCounts {
  sensitive: number;
  resistantA: number;
  resistantB: number;
  resistantBoth: number;
  total: number;
}

export function createCloneState(random: () => number = Math.random): CloneState {
  const cells = new Uint8Array(GRID_WIDTH * GRID_HEIGHT);
  const centerX = GRID_WIDTH / 2;
  const centerY = GRID_HEIGHT / 2;
  for (let y = 0; y < GRID_HEIGHT; y++) {
    for (let x = 0; x < GRID_WIDTH; x++) {
      if (Math.hypot(x - centerX, y - centerY) < 7) {
        cells[y * GRID_WIDTH + x] = 1;
      }
    }
  }
  // A few resistant cells are already there at diagnosis, hidden among the sensitive ones
  const occupied = cells.reduce<number[]>((list, value, index) => (value ? [...list, index] : list), []);
  for (const resistance of [1, 2, 1, 2]) {
    cells[occupied[Math.floor(random() * occupied.length)]] = 1 + resistance;
  }
  return { cells, week: 0, toxicity: 0, outcome: 'playing' };
}

export function countClones(cells: Uint8Array): CloneCounts {
  const counts: CloneCounts = { sensitive: 0, resistantA: 0, resistantB: 0, resistantBoth: 0, total: 0 };
  for (const value of cells) {
    if (value === 1) counts.sensitive++;
    else if (value === 2) counts.resistantA++;
    else if (value === 3) counts.resistantB++;
    else if (value === 4) counts.resistantBoth++;
  }
  counts.total = counts.sensitive + counts.resistantA + counts.resistantB + counts.resistantBoth;
  return counts;
}

const NEIGHBOURS = [[-1, 0], [1, 0], [0, -1], [0, 1], [-1, -1], [1, 1], [-1, 1], [1, -1]];

// One week: drugs kill the cells they can, the survivors divide into free neighbouring space, toxicity changes
export function stepCloneState(state: CloneState, treatment: Treatment, random: () => number = Math.random): CloneState {
  if (state.outcome !== 'playing') {
    return state;
  }
  const cells = state.cells.slice();
  const drugA = treatment === 'A' || treatment === 'AB';
  const drugB = treatment === 'B' || treatment === 'AB';

  for (let index = 0; index < cells.length; index++) {
    const value = cells[index];
    if (!value) continue;
    const bits = value - 1;
    if ((drugA && !(bits & 1) && random() < KILL_RATE) || (drugB && !(bits & 2) && random() < KILL_RATE)) {
      cells[index] = 0;
    }
  }

  // Visit cells in random order so growth has no directional bias
  const order: number[] = [];
  for (let index = 0; index < cells.length; index++) {
    if (cells[index]) order.push(index);
  }
  for (let i = order.length - 1; i > 0; i--) {
    const j = Math.floor(random() * (i + 1));
    [order[i], order[j]] = [order[j], order[i]];
  }
  for (const index of order) {
    const value = cells[index];
    if (!value || random() >= DIVISION_RATE[value]) continue;
    const x = index % GRID_WIDTH;
    const y = Math.floor(index / GRID_WIDTH);
    const free: number[] = [];
    for (const [dx, dy] of NEIGHBOURS) {
      const nx = x + dx;
      const ny = y + dy;
      if (nx >= 0 && nx < GRID_WIDTH && ny >= 0 && ny < GRID_HEIGHT && !cells[ny * GRID_WIDTH + nx]) {
        free.push(ny * GRID_WIDTH + nx);
      }
    }
    if (!free.length) continue;
    let child = value;
    if (random() < MUTATION_RATE) {
      child = 1 + ((value - 1) | (random() < 0.5 ? 1 : 2));
    }
    cells[free[Math.floor(random() * free.length)]] = child;
  }

  const toxicity = Math.max(0, Math.min(100, state.toxicity + TOXICITY_PER_WEEK[treatment]));
  const week = state.week + 1;
  const total = countClones(cells).total;
  let outcome: CloneState['outcome'] = 'playing';
  if (total === 0) outcome = 'cleared';
  else if (total / cells.length >= MAX_BURDEN) outcome = 'progression';
  else if (toxicity >= 100) outcome = 'toxicity';
  else if (week >= WIN_WEEKS) outcome = 'survived';
  return { cells, week, toxicity, outcome };
}
