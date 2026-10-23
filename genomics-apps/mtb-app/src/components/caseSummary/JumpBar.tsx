import { useEffect, useRef, useState } from 'react';
import { ACCENT_CLASSES, CASE_SECTION_ORDER, CASE_SECTIONS, CaseSectionKey, sectionAnchor } from './sections';

// Space between the jump bar and a section scrolled to it
const SCROLL_GAP = 12;
// After a click, keep that section highlighted while the smooth scroll runs
const CLICK_HOLD_MS = 900;

/**
 * A row of links to each Case summary section, pinned just under the patient banner. Clicking one scrolls
 * to that card; the section currently at the top of the view is highlighted.
 */
export default function JumpBar() {
  const barRef = useRef<HTMLElement | null>(null);
  const [bannerHeight, setBannerHeight] = useState(0);
  const [active, setActive] = useState<CaseSectionKey>(CASE_SECTION_ORDER[0]);
  const holdUntilRef = useRef(0);

  // Sit just below the sticky banner, and tell the cards how far below the top to stop when scrolled to
  useEffect(() => {
    const banner = document.getElementById('patient-banner');
    const bar = barRef.current;
    if (!banner || !bar) {
      return;
    }
    const update = () => {
      setBannerHeight(banner.offsetHeight);
      document.documentElement.style.setProperty('--case-scroll-offset', `${banner.offsetHeight + bar.offsetHeight + SCROLL_GAP}px`);
    };
    update();
    const observer = new ResizeObserver(update);
    observer.observe(banner);
    observer.observe(bar);
    return () => {
      observer.disconnect();
      document.documentElement.style.removeProperty('--case-scroll-offset');
    };
  }, []);

  // Highlight the section whose top most recently passed under the jump bar. Cards side by side in a grid
  // row share a top; the first of them wins
  useEffect(() => {
    let frame = 0;
    const onScroll = () => {
      cancelAnimationFrame(frame);
      frame = requestAnimationFrame(() => {
        if (Date.now() < holdUntilRef.current) {
          return;
        }
        const line = (barRef.current?.getBoundingClientRect().bottom ?? 0) + SCROLL_GAP + 4;
        let current = CASE_SECTION_ORDER[0];
        let currentTop = -Infinity;
        for (const key of CASE_SECTION_ORDER) {
          const top = document.getElementById(sectionAnchor(key))?.getBoundingClientRect().top;
          if (top !== undefined && top <= line && top > currentTop + 1) {
            current = key;
            currentTop = top;
          }
        }
        // At the very bottom of the page, the last sections can't reach the line; highlight the last one
        if (window.innerHeight + window.scrollY >= document.documentElement.scrollHeight - 2) {
          current = CASE_SECTION_ORDER[CASE_SECTION_ORDER.length - 1];
        }
        setActive(current);
      });
    };
    onScroll();
    window.addEventListener('scroll', onScroll, { passive: true });
    window.addEventListener('resize', onScroll);
    return () => {
      cancelAnimationFrame(frame);
      window.removeEventListener('scroll', onScroll);
      window.removeEventListener('resize', onScroll);
    };
  }, []);

  const jumpTo = (key: CaseSectionKey) => {
    holdUntilRef.current = Date.now() + CLICK_HOLD_MS;
    document.getElementById(sectionAnchor(key))?.scrollIntoView({ behavior: 'smooth', block: 'start' });
    setActive(key);
  };

  return (
    <nav
      ref={barRef}
      aria-label="Case Summary sections"
      className="sticky z-20 flex flex-wrap gap-1.5 rounded-2xl bg-white/95 p-2 shadow-sm ring-1 ring-inset ring-slate-200 backdrop-blur"
      style={{ top: bannerHeight }}
    >
      {CASE_SECTION_ORDER.map((key) => {
        const { short, icon: Icon, accent } = CASE_SECTIONS[key];
        const isActive = key === active;
        return (
          <button
            key={key}
            type="button"
            onClick={() => jumpTo(key)}
            aria-current={isActive ? 'location' : undefined}
            className={`inline-flex items-center gap-2 rounded-xl px-2.5 py-1.5 text-sm transition-colors ${isActive ? 'bg-slate-100 font-semibold text-slate-900 ring-1 ring-inset ring-slate-300' : 'font-medium text-slate-600 hover:bg-slate-50 hover:text-slate-900'}`}
          >
            <span className={`flex h-6 w-6 items-center justify-center rounded-md ${ACCENT_CLASSES[accent]}`} aria-hidden="true">
              <Icon size={14} strokeWidth={2.25} />
            </span>
            {short}
          </button>
        );
      })}
    </nav>
  );
}
