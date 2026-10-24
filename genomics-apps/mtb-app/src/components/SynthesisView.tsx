import { useEffect, useRef, useState } from 'react';
import ReactMarkdown, { Components } from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { Sparkles } from 'lucide-react';
import { fetchSynthesis } from '@/services/synthesisService';

// The demo's syntheses are pre-generated (mtb-clinical-data-design.md 6.2); these steps stand in for the
// work a live model would do, shown briefly before the stored result appears
const PROGRESS_STEPS = [
  'Reviewing the clinical record…',
  'Reviewing genomic findings…',
  'Checking treatment evidence and clinical trials…',
  'Writing the synthesis…',
];
const STEP_MS = 700;

type SynthesisState =
  | { status: 'idle' }
  | { status: 'generating'; step: number }
  | { status: 'ready'; markdown: string }
  | { status: 'unavailable' }
  | { status: 'error'; message: string };

// Markdown elements styled to match the Case Summary cards
const MARKDOWN_COMPONENTS: Components = {
  h1: ({ children }) => <h1 className="mb-4 text-2xl font-bold tracking-tight text-slate-900">{children}</h1>,
  h2: ({ children }) => <h2 className="mb-3 mt-8 border-b border-slate-200 pb-1.5 text-lg font-semibold text-slate-900">{children}</h2>,
  h3: ({ children }) => <h3 className="mb-2 mt-5 text-base font-semibold text-slate-900">{children}</h3>,
  p: ({ children }) => <p className="my-2 leading-relaxed text-slate-700">{children}</p>,
  ul: ({ children }) => <ul className="my-2 list-disc space-y-1 pl-6 text-slate-700">{children}</ul>,
  ol: ({ children }) => <ol className="my-2 list-decimal space-y-1 pl-6 text-slate-700">{children}</ol>,
  li: ({ children }) => <li className="leading-relaxed">{children}</li>,
  strong: ({ children }) => <strong className="font-semibold text-slate-900">{children}</strong>,
  a: ({ href, children }) => (
    <a href={href} target="_blank" rel="noopener noreferrer" className="font-medium text-blue-600 hover:underline">{children}</a>
  ),
  blockquote: ({ children }) => (
    <blockquote className="my-4 rounded-lg bg-amber-50 px-4 py-1 text-amber-900 ring-1 ring-inset ring-amber-200 [&_p]:text-amber-900">{children}</blockquote>
  ),
  table: ({ children }) => (
    <div className="my-3 overflow-x-auto rounded-lg border border-slate-200">
      <table className="w-full border-collapse text-sm">{children}</table>
    </div>
  ),
  thead: ({ children }) => <thead className="bg-slate-50 text-left text-xs font-semibold uppercase tracking-wide text-slate-500">{children}</thead>,
  tr: ({ children }) => <tr className="border-t border-slate-200 even:bg-slate-50/60">{children}</tr>,
  th: ({ children }) => <th className="px-3 py-2">{children}</th>,
  td: ({ children }) => <td className="px-3 py-2 align-top text-slate-700">{children}</td>,
  hr: () => <hr className="my-6 border-slate-200" />,
  code: ({ children }) => <code className="rounded bg-slate-100 px-1 py-0.5 text-[0.9em]">{children}</code>,
};

export default function SynthesisView({ subjectId }: { subjectId: string }) {
  const [state, setState] = useState<SynthesisState>({ status: 'idle' });
  // Ignores a generation that finishes after the patient has changed
  const runRef = useRef(0);

  useEffect(() => {
    runRef.current += 1;
    setState({ status: 'idle' });
  }, [subjectId]);

  const generate = async () => {
    const run = ++runRef.current;
    setState({ status: 'generating', step: 0 });
    const result = fetchSynthesis(subjectId.trim()).then(
      (markdown) => (markdown === null ? { status: 'unavailable' as const } : { status: 'ready' as const, markdown }),
      (error) => ({ status: 'error' as const, message: error instanceof Error ? error.message : String(error) }),
    );
    for (let step = 1; step < PROGRESS_STEPS.length; step++) {
      await new Promise((resolve) => window.setTimeout(resolve, STEP_MS));
      if (run !== runRef.current) {
        return;
      }
      setState({ status: 'generating', step });
    }
    await new Promise((resolve) => window.setTimeout(resolve, STEP_MS));
    const outcome = await result;
    if (run === runRef.current) {
      setState(outcome);
    }
  };

  const card = 'rounded-2xl border border-slate-200 bg-white p-6 shadow-[0_1px_3px_rgba(15,23,42,0.06),0_4px_12px_rgba(15,23,42,0.04)]';

  if (!subjectId.trim()) {
    return <div className={`${card} text-center text-slate-600`}>Choose a patient in the banner to generate a synthesis.</div>;
  }

  if (state.status === 'ready') {
    return (
      <article className={`${card} mx-auto max-w-5xl text-sm sm:p-8`}>
        <ReactMarkdown remarkPlugins={[remarkGfm]} components={MARKDOWN_COMPONENTS}>{state.markdown}</ReactMarkdown>
      </article>
    );
  }

  return (
    <div className={`${card} mx-auto flex max-w-3xl flex-col items-center gap-5 py-12 text-center`}>
      <span className="flex h-12 w-12 items-center justify-center rounded-xl bg-violet-100 text-violet-700" aria-hidden="true">
        <Sparkles size={24} strokeWidth={2} />
      </span>
      <div>
        <h2 className="text-lg font-semibold text-slate-900">Case synthesis</h2>
        <p className="mt-1 max-w-xl text-sm text-slate-600">
          A summary of the case drawing on the clinical record and genomic findings, with information gaps and treatment options, including recruiting clinical trials.
        </p>
      </div>

      {state.status === 'idle' && (
        <button
          type="button"
          onClick={generate}
          className="inline-flex items-center gap-2 rounded-lg bg-slate-900 px-5 py-2.5 text-sm font-semibold text-white shadow-sm hover:bg-slate-800 focus:outline-none focus:ring-2 focus:ring-slate-400 focus:ring-offset-2"
        >
          <Sparkles size={16} aria-hidden="true" />
          Generate synthesis
        </button>
      )}

      {state.status === 'generating' && (
        <div className="flex w-full max-w-sm flex-col gap-3" role="status" aria-live="polite">
          <div className="h-1.5 overflow-hidden rounded-full bg-slate-100">
            <div
              className="h-full rounded-full bg-violet-500 transition-[width] duration-700 ease-out"
              style={{ width: `${((state.step + 1) / PROGRESS_STEPS.length) * 100}%` }}
            />
          </div>
          <span className="text-sm text-slate-600">{PROGRESS_STEPS[state.step]}</span>
        </div>
      )}

      {state.status === 'unavailable' && (
        <p className="text-sm text-slate-600">Synthesis not available for this patient in the demo.</p>
      )}

      {state.status === 'error' && (
        <div className="flex flex-col items-center gap-3">
          <p className="text-sm text-red-700" title={state.message}>The synthesis could not be loaded.</p>
          <button type="button" onClick={generate} className="rounded-md border border-slate-300 px-3 py-1.5 text-sm font-medium hover:bg-slate-50">Retry</button>
        </div>
      )}

      <p className="text-xs text-slate-500">Demo: synthesis pre-generated for this case.</p>
    </div>
  );
}
