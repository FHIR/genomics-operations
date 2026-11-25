// Each game's best score, kept in this browser only; storage may be unavailable, so failures are ignored

export function readBest(game: string): number {
  try {
    return Number(window.localStorage.getItem(`mtb-games-best-${game}`)) || 0;
  } catch {
    return 0;
  }
}

export function saveBest(game: string, value: number) {
  try {
    window.localStorage.setItem(`mtb-games-best-${game}`, String(value));
  } catch {
    // Not saved; the score still shows for this session
  }
}
