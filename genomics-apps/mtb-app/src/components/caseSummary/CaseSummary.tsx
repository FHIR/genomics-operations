import { ReactNode, useMemo } from 'react';
import { Dna } from 'lucide-react';
import { AllergyIntolerance, Condition, FamilyMemberHistory, FhirResource, Observation } from '@/types/fhir';
import { CrosswalkRow } from '@/cancerFilter/tumorTypeCrosswalk';
import { ClinicalState } from '@/lib/useClinicalRecord';
import { conceptLabel, conditionOnset, effectiveDate, firstCode, formatDate, getMetastasesAfter, hasCode, SNOMED, stageValueLabel } from '@/lib/clinicalSummary';
import {
  bodySiteText,
  CaseSummaryData,
  classifyClinicalRecord,
  codingsText,
  ComorbidityNote,
  histologyText,
  ObservationGroup,
  observationValueText,
  StageGroup,
} from '@/lib/caseSummary';
import { buildTreatmentTimeline } from '@/lib/treatmentTimeline';
import JumpBar from './JumpBar';
import LabsSection from './LabsSection';
import TreatmentTimeline from './TreatmentTimeline';
import ReportsSection from './ReportsSection';
import { Chip, SectionCard } from './SectionCard';

const NO_KNOWN_ALLERGY_CODES = [`${SNOMED}|716186003`, `${SNOMED}|409137002`];
// Earlier values listed after the latest one in the compact sections
const EARLIER_VALUES_SHOWN = 3;

function conditionDate(condition: Condition) {
  return formatDate(condition.onsetDateTime ?? condition.onsetPeriod?.start ?? condition.recordedDate);
}

function Muted({ children }: { children: ReactNode }) {
  return <span className="text-slate-500">{children}</span>;
}

// A list row: name on the left, date (or other detail) right-aligned and muted, then any details below
function Row({ name, aside, children }: { name: ReactNode; aside?: ReactNode; children?: ReactNode }) {
  return (
    <li className="border-b border-slate-100 py-2 first:pt-0 last:border-b-0">
      <div className="flex items-baseline justify-between gap-3">
        <span className="font-semibold text-slate-900">{name}</span>
        {aside && <span className="shrink-0 text-xs text-slate-500">{aside}</span>}
      </div>
      {children}
    </li>
  );
}

function StageLabel({ observation }: { observation: Observation }) {
  const value = conceptLabel(observation.valueCodeableConcept);
  return <>{value ? `Stage ${stageValueLabel(value)}` : 'Stage not recorded'}</>;
}

function StageCard({ group, laterMetastases = [] }: { group: StageGroup; laterMetastases?: Condition[] }) {
  return (
    <div className="flex flex-col gap-2 rounded-xl bg-slate-50 p-4 text-sm ring-1 ring-inset ring-slate-200">
      <div className="flex items-baseline justify-between gap-3">
        <span className="text-base font-bold text-slate-900"><StageLabel observation={group.observation} /></span>
        <span className="text-xs text-slate-500">{formatDate(effectiveDate(group.observation)) ?? 'no date'}</span>
      </div>
      <Muted>{[conceptLabel(group.observation.code), conceptLabel(group.observation.method)].filter(Boolean).join(' · ')}</Muted>
      {group.members.length > 0 && (
        <div className="flex flex-wrap gap-1.5">
          {group.members.map((member) => {
            const memberValue = conceptLabel(member.valueCodeableConcept);
            return <Chip key={member.id} title={conceptLabel(member.code)}>{memberValue ? stageValueLabel(memberValue) : conceptLabel(member.code)}</Chip>;
          })}
        </div>
      )}
      {laterMetastases.length > 0 && (
        <span className="rounded-lg bg-amber-50 px-2.5 py-1.5 text-amber-900 ring-1 ring-inset ring-amber-200">
          Metastases recorded after this staging: {laterMetastases
            .map((metastasis) => `${conceptLabel(metastasis.code) ?? 'Secondary cancer'} (${formatDate(conditionOnset(metastasis))})`)
            .join('; ')}
        </span>
      )}
    </div>
  );
}

// The feature block: the cancer as recorded, its newest stage as a dark pill, and its details as chips
function PrimaryCancerHeader({ condition, latestStage }: { condition: Condition; latestStage?: Observation }) {
  const histology = histologyText(condition);
  const sites = bodySiteText(condition);
  const clinicalStatus = firstCode(condition.clinicalStatus);
  const verificationStatus = firstCode(condition.verificationStatus);
  return (
    <div className="flex flex-col gap-3 rounded-xl bg-gradient-to-br from-slate-50 to-white p-4 ring-1 ring-inset ring-slate-200">
      <div className="flex flex-wrap items-center gap-3">
        <span className="text-xl font-bold tracking-tight text-slate-900">{conceptLabel(condition.code) ?? 'Primary cancer (no code)'}</span>
        {latestStage && <Chip tone="dark"><StageLabel observation={latestStage} /></Chip>}
      </div>
      <div className="flex flex-wrap gap-1.5">
        {histology && <Chip title="Histology">{histology}</Chip>}
        {sites.map((site) => <Chip key={site} title="Body site">{site}</Chip>)}
        {conditionDate(condition) && <Chip tone="muted">Onset {conditionDate(condition)}</Chip>}
        {clinicalStatus && <Chip tone={['recurrence', 'relapse'].includes(clinicalStatus) ? 'high' : 'muted'} title="Clinical status">{clinicalStatus}</Chip>}
        {verificationStatus && <Chip tone="muted" title="Verification status">{verificationStatus}</Chip>}
      </div>
      <span className="text-xs text-slate-500">{codingsText(condition.code)}</span>
    </div>
  );
}

function DiagnosisSection({ data }: { data: CaseSummaryData }) {
  const hasAnything = data.primaryCancers.length > 0 || data.otherStaging.length > 0 || data.looseTnmCategories.length > 0 || data.metastases.length > 0;
  return (
    <SectionCard section="diagnosis">
      {!hasAnything ? <p className="text-sm text-slate-500">No cancer condition recognized in this record (see Other conditions).</p> : (
        <div className="flex flex-col gap-4">
          {data.primaryCancers.map(({ condition, stageGroups }) => (
            <div key={condition.id} className="flex flex-col gap-3">
              <PrimaryCancerHeader condition={condition} latestStage={stageGroups[0]?.observation} />
              <div className="grid grid-cols-1 gap-3 md:grid-cols-2 xl:grid-cols-3">
                {/* Only the newest staging gets the note; AJCC stage isn't reassigned when disease spreads */}
                {stageGroups.map((group, position) => (
                  <StageCard
                    key={group.observation.id}
                    group={group}
                    laterMetastases={position === 0 ? getMetastasesAfter(data.metastases, effectiveDate(group.observation)) : []}
                  />
                ))}
                {stageGroups.length === 0 && condition.stage?.some((stage) => stage.summary) && (
                  <div className="flex flex-col gap-1 rounded-xl bg-slate-50 p-4 text-sm ring-1 ring-inset ring-slate-200">
                    <span className="text-xs font-semibold uppercase tracking-wide text-slate-500">Stage (recorded on the condition)</span>
                    {condition.stage.map((stage, index) => {
                      const summary = conceptLabel(stage.summary);
                      return summary && <span key={index} className="text-base font-bold text-slate-900">Stage {stageValueLabel(summary)}</span>;
                    })}
                  </div>
                )}
              </div>
            </div>
          ))}
          {(data.otherStaging.length > 0 || data.looseTnmCategories.length > 0) && (
            <div className="grid grid-cols-1 gap-3 md:grid-cols-2 xl:grid-cols-3">
              {data.otherStaging.map((group) => <StageCard key={group.observation.id} group={group} />)}
              {data.looseTnmCategories.length > 0 && (
                <div className="flex flex-col gap-1 rounded-xl bg-slate-50 p-4 text-sm ring-1 ring-inset ring-slate-200">
                  <span className="text-xs font-semibold uppercase tracking-wide text-slate-500">TNM categories not in a stage group</span>
                  {data.looseTnmCategories.map((observation) => (
                    <span key={observation.id}>
                      {conceptLabel(observation.code)}: <b>{stageValueLabel(conceptLabel(observation.valueCodeableConcept) ?? '—')}</b>
                      <Muted> · {formatDate(effectiveDate(observation)) ?? 'no date'}</Muted>
                    </span>
                  ))}
                </div>
              )}
            </div>
          )}
          {data.metastases.length > 0 && (
            <div className="flex flex-col gap-2">
              <span className="text-xs font-semibold uppercase tracking-wide text-slate-500">Metastases</span>
              <div className="flex flex-wrap gap-2">
                {data.metastases.map((condition) => {
                  const sites = bodySiteText(condition);
                  return (
                    <span
                      key={condition.id}
                      title={sites.length > 0 ? `Body site: ${sites.join(', ')}` : undefined}
                      className="inline-flex items-baseline gap-2 rounded-lg bg-white px-3 py-1.5 text-sm ring-1 ring-inset ring-slate-200"
                    >
                      <span className="font-medium text-slate-900">{conceptLabel(condition.code) ?? 'Secondary cancer (no code)'}</span>
                      <span className="text-xs text-slate-500">{conditionDate(condition) ?? 'no date'}</span>
                    </span>
                  );
                })}
              </div>
            </div>
          )}
        </div>
      )}
    </SectionCard>
  );
}

// One line per test: latest value and date, then a few earlier values
function ObservationGroupList({ groups, index }: { groups: ObservationGroup[]; index: Map<string, FhirResource> }) {
  if (groups.length === 0) {
    return <p className="text-sm text-slate-500">None recorded</p>;
  }
  const valueAndDate = (observation: Observation) =>
    `${observationValueText(observation, index) ?? '—'} · ${formatDate(effectiveDate(observation)) ?? 'no date'}`;
  return (
    <ul className="text-sm">
      {groups.map((group) => {
        const [latest, ...earlier] = group.observations;
        return (
          <Row key={group.key} name={group.label} aside={formatDate(effectiveDate(latest)) ?? 'no date'}>
            <span className="text-slate-800">{observationValueText(latest, index) ?? '—'}</span>
            {earlier.length > 0 && (
              <span className="mt-0.5 block text-xs text-slate-500">
                Earlier: {earlier.slice(0, EARLIER_VALUES_SHOWN).map(valueAndDate).join('; ')}
                {earlier.length > EARLIER_VALUES_SHOWN && `; +${earlier.length - EARLIER_VALUES_SHOWN} more`}
              </span>
            )}
          </Row>
        );
      })}
    </ul>
  );
}

interface OtherConditionsListProps {
  conditions: Condition[];
  // From comorbidities assessments, newest first
  absent: ComorbidityNote[];
  notListed: ComorbidityNote[];
}

function OtherConditionsList({ conditions, absent, notListed }: OtherConditionsListProps) {
  if (conditions.length === 0 && absent.length === 0 && notListed.length === 0) {
    return <p className="text-sm text-slate-500">None recorded</p>;
  }
  return (
    <ul className="text-sm">
      {conditions.map((condition) => (
        <Row
          key={condition.id}
          name={conceptLabel(condition.code) ?? 'Condition (no code)'}
          aside={[firstCode(condition.clinicalStatus), conditionDate(condition) && `since ${conditionDate(condition)}`].filter(Boolean).join(' · ')}
        />
      ))}
      {notListed.map((note) => (
        <Row key={`noted-${note.date}`} name="Also noted as comorbidities" aside={note.date ?? 'no date'}>
          <span className="text-slate-700">{note.labels.join(', ')}</span>
        </Row>
      ))}
      {absent.map((note) => (
        <Row key={`absent-${note.date}`} name="Documented absent" aside={note.date ?? 'no date'}>
          <span className="text-slate-700">{note.labels.join(', ')}</span>
        </Row>
      ))}
    </ul>
  );
}

function AllergyList({ allergies }: { allergies: AllergyIntolerance[] }) {
  // No AllergyIntolerance at all isn't the same as a recorded "no known allergies"
  if (allergies.length === 0) {
    return <p className="text-sm text-slate-500">None recorded</p>;
  }
  return (
    <ul className="text-sm">
      {allergies.map((allergy) => {
        if (hasCode(allergy.code, NO_KNOWN_ALLERGY_CODES)) {
          return <Row key={allergy.id} name="No known allergies" />;
        }
        const status = firstCode(allergy.clinicalStatus);
        const isInactive = status === 'inactive' || status === 'resolved';
        const reactions = (allergy.reaction ?? []).map((reaction) => {
          const manifestations = (reaction.manifestation ?? []).map((manifestation) => conceptLabel(manifestation)).filter(Boolean).join(', ');
          return reaction.severity ? `${manifestations || 'Reaction'} (${reaction.severity})` : manifestations;
        }).filter(Boolean);
        // High criticality stands out as a red pill; other details stay plain text
        const isHighCriticality = allergy.criticality === 'high';
        const details = [
          reactions.join('; '),
          !isHighCriticality && allergy.criticality && `criticality ${allergy.criticality}`,
          allergy.category?.join(', '),
        ].filter(Boolean);
        return (
          <Row
            key={allergy.id}
            name={<span className={isInactive ? 'text-slate-500' : ''}>{conceptLabel(allergy.code) ?? 'Unnamed allergy'}</span>}
            aside={status && status !== 'active' ? status : undefined}
          >
            {(details.length > 0 || isHighCriticality) && (
              <div className={`mt-0.5 flex flex-wrap items-center gap-2 ${isInactive ? 'text-slate-500' : 'text-slate-700'}`}>
                {isHighCriticality && <Chip tone="danger">High criticality</Chip>}
                {details.length > 0 && <span>{details.join(' · ')}</span>}
              </div>
            )}
          </Row>
        );
      })}
    </ul>
  );
}

// One line per relative: "Mother: Malignant neoplasm of breast (age 64)"
function FamilyHistoryList({ entries }: { entries: FamilyMemberHistory[] }) {
  if (entries.length === 0) {
    return <p className="text-sm text-slate-500">None recorded</p>;
  }
  return (
    <ul className="text-sm">
      {entries.map((entry) => {
        const relationship = conceptLabel(entry.relationship) ?? 'Relative';
        const conditions = (entry.condition ?? []).map((condition) => {
          const age = condition.onsetAge?.value !== undefined ? `age ${condition.onsetAge.value}` : condition.onsetString;
          const details = [age, condition.contributedToDeath && 'contributed to death'].filter(Boolean).join(', ');
          return `${conceptLabel(condition.code) ?? 'Condition'}${details ? ` (${details})` : ''}`;
        });
        return (
          <Row key={entry.id} name={relationship.charAt(0).toUpperCase() + relationship.slice(1)} aside={entry.deceasedBoolean ? 'deceased' : undefined}>
            {conditions.length > 0 && <span className="text-slate-700">{conditions.join('; ')}</span>}
          </Row>
        );
      })}
    </ul>
  );
}

interface CaseSummaryProps {
  subjectId: string;
  clinicalState: ClinicalState;
  onRetry: () => void;
  // Offered when the patient has no clinical data
  onGoToGenomicFindings: () => void;
  crosswalk: CrosswalkRow[];
}

export default function CaseSummary({ subjectId, clinicalState, onRetry, onGoToGenomicFindings, crosswalk }: CaseSummaryProps) {
  const record = clinicalState.status === 'loaded' ? clinicalState.record : null;
  const data = useMemo(() => (record ? classifyClinicalRecord(record, crosswalk) : null), [record, crosswalk]);
  const timeline = useMemo(() => {
    if (!data) {
      return null;
    }
    // Today in local time, as YYYY-MM-DD
    const now = new Date();
    const today = `${now.getFullYear()}-${String(now.getMonth() + 1).padStart(2, '0')}-${String(now.getDate()).padStart(2, '0')}`;
    return buildTreatmentTimeline(data, today);
  }, [data]);

  if (!data) {
    return (
      <div className="rounded-2xl border border-slate-200 bg-white p-8 text-center text-slate-600 shadow-sm">
        {!subjectId.trim() && 'Choose a patient in the banner to see their clinical data.'}
        {clinicalState.status === 'loading' && 'Loading clinical data…'}
        {clinicalState.status === 'none' && (
          <span className="inline-flex flex-col items-center gap-3">
            <span>No clinical data found for MRN {subjectId}. Genomic Findings still work as usual.</span>
            <button
              type="button"
              onClick={onGoToGenomicFindings}
              className="inline-flex items-center gap-2 rounded-lg bg-slate-900 px-4 py-2 text-sm font-semibold text-white hover:bg-slate-800"
            >
              <Dna size={16} aria-hidden="true" className="text-emerald-300" />
              Go to Genomic Findings
            </button>
          </span>
        )}
        {clinicalState.status === 'error' && (
          <span className="inline-flex flex-col items-center gap-3">
            <span className="text-red-700" title={clinicalState.message}>Clinical data could not be loaded.</span>
            <button type="button" onClick={onRetry} className="rounded-md border border-slate-300 px-3 py-1.5 text-sm font-medium hover:bg-slate-50">Retry</button>
          </span>
        )}
      </div>
    );
  }

  return (
    <div className="flex flex-col gap-5">
      <JumpBar />
      <DiagnosisSection data={data} />
      <TreatmentTimeline timeline={timeline} />
      <LabsSection labs={data.labs} index={data.index} />
      {/* Three cards per row on wide screens, two on medium, one on narrow */}
      <div className="grid grid-cols-1 gap-5 md:grid-cols-2 xl:grid-cols-3">
        <SectionCard section="imaging" scrollBody>
          <ObservationGroupList groups={data.imaging} index={data.index} />
        </SectionCard>
        <SectionCard section="assessments" scrollBody>
          <ObservationGroupList groups={data.assessments} index={data.index} />
        </SectionCard>
        <SectionCard section="vitals" scrollBody>
          <ObservationGroupList groups={data.vitals} index={data.index} />
        </SectionCard>
        <SectionCard section="conditions" scrollBody>
          <OtherConditionsList conditions={data.otherConditions} absent={data.comorbiditiesAbsent} notListed={data.comorbiditiesNotListed} />
        </SectionCard>
        <SectionCard section="allergies" scrollBody>
          <AllergyList allergies={data.allergies} />
        </SectionCard>
        <SectionCard section="family" scrollBody>
          <FamilyHistoryList entries={data.familyHistory} />
        </SectionCard>
      </div>
      <ReportsSection groups={data.reportGroups} />
    </div>
  );
}
