import { useEffect, useMemo, useState } from 'react';
import { CaseDx } from '@/cancerFilter/tumorTypeCrosswalk';
import { fetchClinicalMrns } from '@/services/clinicalService';
import { ClinicalState } from '@/lib/useClinicalRecord';
import {
  DiseaseStatusTrend,
  getDiseaseStatus,
  getLatestEcog,
  getLatestStage,
  getPatientSummary,
} from '@/lib/clinicalSummary';
import { Condition } from '@/types/fhir';

// Patients with genomic data on the genomics operations server
export const GENOMIC_PATIENT_IDS: string[] = [
  "ABC123", "NA19240", "ABC456",
  "NA18498", "NA19247", "m123",
  "ABC789", "NA18499", "NA19256",
  "CA12345", "NA18870", "NB6TK328",
  "HCC1143", "NA18871", "NB6TK329",
  "HG00403", "NA19190", "XYZ123",
  "HG00406", "NA19210", "XYZ234",
  "HG02657", "NA19238", "XYZ345",
  "huC30902", "NA19239", "L2345",
  "MRN1234", "TCGA-DD-A1EH"
];

const CUSTOM_MRN_VALUE = '__custom__';

function sortIds(ids: string[]) {
  return [...ids].sort((a, b) => a.toLowerCase().localeCompare(b.toLowerCase()));
}

// Filled chips that read on the dark banner
const STATUS_CHIP_CLASSES: Record<DiseaseStatusTrend, string> = {
  progression: 'bg-orange-300 text-orange-950',
  improving: 'bg-emerald-300 text-emerald-950',
  stable: 'bg-slate-300 text-slate-900',
  other: 'bg-amber-200 text-amber-950',
};

// "male" -> "M"; other values (e.g. "unknown") are shown as recorded
const SEX_ABBREVIATIONS: Record<string, string> = { male: 'M', female: 'F' };

// One line at desktop width: long Dx labels truncate in the closed dropdown (full label in its tooltip)
const selectClassName = 'max-w-[18rem] truncate rounded-md border border-slate-500 bg-slate-800 px-2 py-1 text-sm font-semibold text-white focus:outline-none focus:ring-2 focus:ring-sky-400';

interface MrnSelectProps {
  subjectId: string;
  onSubjectIdChange: (value: string) => void;
}

function MrnSelect({ subjectId, onSubjectIdChange }: MrnSelectProps) {
  const [isCustomMrn, setIsCustomMrn] = useState(false);
  // null until detection finishes; stays null if it fails, and the list is then shown ungrouped
  const [clinicalMrns, setClinicalMrns] = useState<string[] | null>(null);

  useEffect(() => {
    let isMounted = true;
    fetchClinicalMrns()
      .then((mrns) => {
        if (isMounted) {
          setClinicalMrns(mrns);
        }
      })
      .catch((error) => console.error('Error detecting patients with clinical data:', error));
    return () => {
      isMounted = false;
    };
  }, []);

  const withClinicalData = useMemo(() => sortIds(clinicalMrns ?? []), [clinicalMrns]);
  const genomicOnly = useMemo(
    () => sortIds(GENOMIC_PATIENT_IDS.filter((id) => !withClinicalData.includes(id))),
    [withClinicalData]
  );

  if (isCustomMrn) {
    return (
      <span className="inline-flex items-center gap-1">
        <input
          type="text"
          value={subjectId}
          onChange={(e) => onSubjectIdChange(e.target.value)}
          placeholder="Enter custom MRN"
          aria-label="Custom MRN"
          autoFocus
          className="w-40 rounded-md border border-slate-500 bg-slate-800 px-2 py-1 text-sm font-semibold text-white placeholder:text-slate-400 focus:outline-none focus:ring-2 focus:ring-sky-400"
        />
        <button
          type="button"
          onClick={() => {
            setIsCustomMrn(false);
            onSubjectIdChange('');
          }}
          className="rounded-md border border-slate-500 px-2 py-1 text-xs text-slate-200 hover:bg-slate-700"
        >
          Back to list
        </button>
      </span>
    );
  }

  const renderOptions = (ids: string[]) => ids.map((id) => <option key={id} value={id}>{id}</option>);

  return (
    <select
      aria-label="MRN"
      value={subjectId}
      onChange={(e) => {
        if (e.target.value === CUSTOM_MRN_VALUE) {
          setIsCustomMrn(true);
          onSubjectIdChange('');
        } else {
          onSubjectIdChange(e.target.value);
        }
      }}
      className={selectClassName}
    >
      <option value="">Select a patient</option>
      {clinicalMrns === null ? renderOptions(sortIds(GENOMIC_PATIENT_IDS)) : (
        <>
          <optgroup label="With clinical data">{renderOptions(withClinicalData)}</optgroup>
          <optgroup label="Genomic data only">{renderOptions(genomicOnly)}</optgroup>
        </>
      )}
      <option value={CUSTOM_MRN_VALUE}>Enter custom MRN…</option>
    </select>
  );
}

interface DxSelectProps {
  patientOptions: CaseDx[];
  genericOptions: CaseDx[];
  selectedDx: CaseDx | null;
  onDxChange: (dx: CaseDx | null) => void;
}

function DxSelect({ patientOptions, genericOptions, selectedDx, onDxChange }: DxSelectProps) {
  const options = [...patientOptions, ...genericOptions];
  return (
    <select
      aria-label="Diagnosis"
      value={selectedDx?.key ?? ''}
      onChange={(e) => onDxChange(options.find((option) => option.key === e.target.value) ?? null)}
      title={selectedDx?.label}
      className={selectClassName}
    >
      <option value="">No Dx selected</option>
      {patientOptions.map((option) => <option key={option.key} value={option.key}>{option.label}</option>)}
      {patientOptions.length > 0 && genericOptions.length > 0 && <option disabled>──────────</option>}
      {genericOptions.map((option) => <option key={option.key} value={option.key}>{option.label}</option>)}
    </select>
  );
}

interface PatientBannerProps extends MrnSelectProps, DxSelectProps {
  clinicalState: ClinicalState;
  onRetry: () => void;
  // The primary cancer that stage and disease status are read for
  bannerCondition: Condition | undefined;
}

function Caption({ children }: { children: React.ReactNode }) {
  return <span className="text-slate-400">{children}</span>;
}

export default function PatientBanner({
  subjectId,
  onSubjectIdChange,
  patientOptions,
  genericOptions,
  selectedDx,
  onDxChange,
  clinicalState,
  onRetry,
  bannerCondition,
}: PatientBannerProps) {
  const record = clinicalState.status === 'loaded' ? clinicalState.record : null;
  const summary = record ? getPatientSummary(record.patient) : null;
  const stage = record ? getLatestStage(record.resources, bannerCondition) : undefined;
  const ecog = record ? getLatestEcog(record.resources) : undefined;
  const status = record ? getDiseaseStatus(record.resources, bannerCondition) : undefined;

  const demographics = summary ? [
    summary.age !== undefined ? `${summary.age} y` : undefined,
    summary.sex ? (SEX_ABBREVIATIONS[summary.sex] ?? summary.sex) : undefined,
    summary.birthDate ? `DOB ${summary.birthDate}` : undefined,
  ].filter(Boolean) : [];
  const clinicalItems = [stage, ecog].filter((item) => item !== undefined);

  return (
    <div className="sticky top-0 z-30 -ml-64 mb-6 flex w-full min-w-[1400px] flex-nowrap items-center gap-x-4 whitespace-nowrap rounded-xl bg-slate-900 px-6 py-3 text-sm text-slate-200 shadow-md">
      {summary?.name && <span className="text-lg font-bold text-white">{summary.name}</span>}
      {demographics.length > 0 && <span className="text-slate-300">{demographics.join(' · ')}</span>}

      <label className="inline-flex items-center gap-2">
        <Caption>MRN</Caption>
        <MrnSelect subjectId={subjectId} onSubjectIdChange={onSubjectIdChange} />
      </label>

      {subjectId.trim() && (
        <label className="inline-flex items-center gap-2">
          <Caption>Dx</Caption>
          <DxSelect
            patientOptions={patientOptions}
            genericOptions={genericOptions}
            selectedDx={selectedDx}
            onDxChange={onDxChange}
          />
        </label>
      )}

      {clinicalItems.map((item) => (
        <span key={item.caption} title={item.tooltip} className="inline-flex items-center gap-1.5">
          <Caption>{item.caption}</Caption>
          <span className="font-bold text-white">{item.value}</span>
          {item.note && <span className="text-slate-400">· {item.note}</span>}
        </span>
      ))}

      {status && (
        <span
          title={status.tooltip}
          className={`shrink-0 rounded-full px-3 py-1 text-xs font-bold ${STATUS_CHIP_CLASSES[status.trend]}`}
        >
          {status.label}{status.date ? ` ${status.date}` : ''}
        </span>
      )}

      {!subjectId.trim() && <span className="text-slate-400">Choose a patient to begin.</span>}
      {clinicalState.status === 'loading' && <span className="text-slate-400">Loading clinical data…</span>}
      {clinicalState.status === 'none' && <span className="text-slate-400">No clinical data found</span>}
      {clinicalState.status === 'error' && (
        <span className="inline-flex items-center gap-2 text-red-300" title={clinicalState.message}>
          Clinical data unavailable
          <button
            type="button"
            onClick={onRetry}
            className="rounded-md border border-red-400 px-2 py-0.5 text-xs font-medium hover:bg-red-950"
          >
            Retry
          </button>
        </span>
      )}
    </div>
  );
}
