import { useMemo, useState } from 'react';
import { AllergyIntolerance, Condition, FhirResource, Observation } from '@/types/fhir';
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
  resourceLabel,
  StageGroup,
} from '@/lib/caseSummary';
import LabsSection from './LabsSection';
import ReportsSection from './ReportsSection';
import { SectionCard } from './SectionCard';

const NO_KNOWN_ALLERGY_CODES = [`${SNOMED}|716186003`, `${SNOMED}|409137002`];
// Earlier values listed after the latest one in the compact sections
const EARLIER_VALUES_SHOWN = 3;

function conditionDate(condition: Condition) {
  return formatDate(condition.onsetDateTime ?? condition.onsetPeriod?.start ?? condition.recordedDate);
}

function Muted({ children }: { children: React.ReactNode }) {
  return <span className="text-slate-500">{children}</span>;
}

function StageCard({ group, laterMetastases = [] }: { group: StageGroup; laterMetastases?: Condition[] }) {
  const value = conceptLabel(group.observation.valueCodeableConcept);
  return (
    <div className="flex flex-col gap-0.5 rounded-lg border border-slate-200 p-3 text-sm">
      <span className="text-xs font-semibold uppercase tracking-wide text-slate-500">
        Stage · {formatDate(effectiveDate(group.observation)) ?? 'no date'}
      </span>
      <span className="text-base font-semibold text-slate-900">{value ? `Stage ${stageValueLabel(value)}` : 'Stage not recorded'}</span>
      <Muted>{[conceptLabel(group.observation.code), conceptLabel(group.observation.method)].filter(Boolean).join(' · ')}</Muted>
      {group.members.length > 0 && (
        <span className="text-slate-700">
          {group.members.map((member) => {
            const memberValue = conceptLabel(member.valueCodeableConcept);
            return memberValue ? stageValueLabel(memberValue) : conceptLabel(member.code);
          }).join(' · ')}
        </span>
      )}
      {laterMetastases.length > 0 && (
        <span className="mt-1 rounded-md bg-amber-50 px-2 py-1 text-amber-900">
          Metastases recorded after this staging: {laterMetastases
            .map((metastasis) => `${conceptLabel(metastasis.code) ?? 'Secondary cancer'} (${formatDate(conditionOnset(metastasis))})`)
            .join('; ')}
        </span>
      )}
    </div>
  );
}

function PrimaryCancerCard({ condition }: { condition: Condition }) {
  const histology = histologyText(condition);
  const sites = bodySiteText(condition);
  const status = [conditionDate(condition) && `Onset ${conditionDate(condition)}`, firstCode(condition.clinicalStatus), firstCode(condition.verificationStatus)]
    .filter(Boolean).join(' · ');
  return (
    <div className="flex flex-col gap-0.5 rounded-lg border border-slate-200 p-3 text-sm">
      <span className="text-xs font-semibold uppercase tracking-wide text-slate-500">Primary cancer</span>
      <span className="text-base font-semibold text-slate-900">{conceptLabel(condition.code) ?? 'Primary cancer (no code)'}</span>
      <Muted>{codingsText(condition.code)}</Muted>
      {histology && <span className="text-slate-700">Histology: {histology}</span>}
      {sites.length > 0 && <span className="text-slate-700">Body site: {sites.join(', ')}</span>}
      {status && <span className="text-slate-700">{status}</span>}
    </div>
  );
}

function DiagnosisSection({ data }: { data: CaseSummaryData }) {
  const hasAnything = data.primaryCancers.length > 0 || data.otherStaging.length > 0 || data.looseTnmCategories.length > 0 || data.metastases.length > 0;
  return (
    <SectionCard title="Cancer diagnosis & staging">
      {!hasAnything ? <p className="text-sm text-slate-500">No cancer condition recognized in this record (see Other conditions).</p> : (
        <div className="flex flex-col gap-4">
          {data.primaryCancers.map(({ condition, stageGroups }) => (
            <div key={condition.id} className="grid grid-cols-1 gap-3 lg:grid-cols-3">
              <PrimaryCancerCard condition={condition} />
              {/* Only the newest staging gets the note; AJCC stage isn't reassigned when disease spreads */}
              {stageGroups.map((group, position) => (
                <StageCard
                  key={group.observation.id}
                  group={group}
                  laterMetastases={position === 0 ? getMetastasesAfter(data.metastases, effectiveDate(group.observation)) : []}
                />
              ))}
              {stageGroups.length === 0 && condition.stage?.some((stage) => stage.summary) && (
                <div className="flex flex-col gap-0.5 rounded-lg border border-slate-200 p-3 text-sm">
                  <span className="text-xs font-semibold uppercase tracking-wide text-slate-500">Stage (recorded on the condition)</span>
                  {condition.stage.map((stage, index) => {
                    const summary = conceptLabel(stage.summary);
                    return summary && <span key={index} className="text-base font-semibold text-slate-900">Stage {stageValueLabel(summary)}</span>;
                  })}
                </div>
              )}
            </div>
          ))}
          {(data.otherStaging.length > 0 || data.looseTnmCategories.length > 0) && (
            <div className="grid grid-cols-1 gap-3 lg:grid-cols-3">
              {data.otherStaging.map((group) => <StageCard key={group.observation.id} group={group} />)}
              {data.looseTnmCategories.length > 0 && (
                <div className="flex flex-col gap-0.5 rounded-lg border border-slate-200 p-3 text-sm">
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
            <div className="flex flex-col gap-1 rounded-lg border border-slate-200 p-3 text-sm">
              <span className="text-xs font-semibold uppercase tracking-wide text-slate-500">Metastases</span>
              {data.metastases.map((condition) => {
                const sites = bodySiteText(condition);
                return (
                  <span key={condition.id}>
                    <b className="text-slate-900">{conceptLabel(condition.code) ?? 'Secondary cancer (no code)'}</b>
                    {sites.length > 0 && <span className="text-slate-700"> · {sites.join(', ')}</span>}
                    <Muted> · {conditionDate(condition) ?? 'no date'}</Muted>
                  </span>
                );
              })}
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
    <ul className="flex flex-col gap-2 text-sm">
      {groups.map((group) => {
        const [latest, ...earlier] = group.observations;
        return (
          <li key={group.key}>
            <b className="text-slate-900">{group.label}</b>
            <br />
            <span className="text-slate-800">{valueAndDate(latest)}</span>
            {earlier.length > 0 && (
              <span className="block text-xs text-slate-500">
                Earlier: {earlier.slice(0, EARLIER_VALUES_SHOWN).map(valueAndDate).join('; ')}
                {earlier.length > EARLIER_VALUES_SHOWN && `; +${earlier.length - EARLIER_VALUES_SHOWN} more`}
              </span>
            )}
          </li>
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
    <ul className="flex flex-col gap-1 text-sm">
      {conditions.map((condition) => (
        <li key={condition.id}>
          <b className="text-slate-900">{conceptLabel(condition.code) ?? 'Condition (no code)'}</b>
          <Muted>
            {[firstCode(condition.clinicalStatus), conditionDate(condition) && `since ${conditionDate(condition)}`].filter(Boolean).map((part) => ` · ${part}`).join('')}
          </Muted>
        </li>
      ))}
      {notListed.map((note) => (
        <li key={`noted-${note.date}`}><span className="text-slate-700">Also noted as comorbidities</span><Muted> ({note.date ?? 'no date'})</Muted>: {note.labels.join(', ')}</li>
      ))}
      {absent.map((note) => (
        <li key={`absent-${note.date}`}><span className="text-slate-700">Documented absent</span><Muted> ({note.date ?? 'no date'})</Muted>: {note.labels.join(', ')}</li>
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
    <ul className="flex flex-col gap-1 text-sm">
      {allergies.map((allergy) => {
        if (hasCode(allergy.code, NO_KNOWN_ALLERGY_CODES)) {
          return <li key={allergy.id}><b className="text-slate-900">No known allergies</b></li>;
        }
        const status = firstCode(allergy.clinicalStatus);
        const isInactive = status === 'inactive' || status === 'resolved';
        const reactions = (allergy.reaction ?? []).map((reaction) => {
          const manifestations = (reaction.manifestation ?? []).map((manifestation) => conceptLabel(manifestation)).filter(Boolean).join(', ');
          return reaction.severity ? `${manifestations || 'Reaction'} (${reaction.severity})` : manifestations;
        }).filter(Boolean);
        const details = [
          reactions.join('; '),
          allergy.criticality && `criticality ${allergy.criticality}`,
          allergy.category?.join(', '),
          status && status !== 'active' && status,
        ].filter(Boolean);
        return (
          <li key={allergy.id} className={isInactive ? 'text-slate-500' : ''}>
            <b className={isInactive ? 'text-slate-600' : 'text-slate-900'}>{conceptLabel(allergy.code) ?? 'Unnamed allergy'}</b>
            {details.length > 0 && <span className={isInactive ? '' : 'text-slate-700'}> · {details.join(' · ')}</span>}
          </li>
        );
      })}
    </ul>
  );
}

function OtherResources({ groups }: { groups: CaseSummaryData['otherResources'] }) {
  const [openType, setOpenType] = useState<string | null>(null);
  const open = groups.find((group) => group.resourceType === openType);
  return (
    <SectionCard title="Other resources">
      {groups.length === 0 ? <p className="text-sm text-slate-500">None</p> : (
        <>
          <div className="flex flex-wrap gap-2">
            {groups.map((group) => (
              <button
                key={group.resourceType}
                type="button"
                aria-expanded={openType === group.resourceType}
                onClick={() => setOpenType(openType === group.resourceType ? null : group.resourceType)}
                className={`rounded-full border px-3 py-1 text-sm font-medium ${openType === group.resourceType ? 'border-blue-700 bg-blue-600 text-white' : 'border-slate-300 bg-white text-slate-700 hover:bg-slate-50'}`}
              >
                {group.resourceType} ({group.resources.length})
              </button>
            ))}
          </div>
          {open && (
            <ul className="mt-3 flex flex-col gap-1 border-t border-slate-200 pt-3 text-sm">
              {open.resources.map((resource) => (
                <li key={resource.id}>{resourceLabel(resource) || resource.resourceType}</li>
              ))}
            </ul>
          )}
        </>
      )}
    </SectionCard>
  );
}

interface CaseSummaryProps {
  subjectId: string;
  clinicalState: ClinicalState;
  onRetry: () => void;
  crosswalk: CrosswalkRow[];
}

export default function CaseSummary({ subjectId, clinicalState, onRetry, crosswalk }: CaseSummaryProps) {
  const record = clinicalState.status === 'loaded' ? clinicalState.record : null;
  const data = useMemo(() => (record ? classifyClinicalRecord(record, crosswalk) : null), [record, crosswalk]);

  if (!data) {
    return (
      <div className="rounded-xl border border-slate-300 bg-white p-8 text-center text-slate-600 shadow-sm">
        {!subjectId.trim() && 'Choose a patient in the banner to see their clinical data.'}
        {clinicalState.status === 'loading' && 'Loading clinical data…'}
        {clinicalState.status === 'none' && `No clinical data found for MRN ${subjectId}. Genomic findings still work as usual.`}
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
      <DiagnosisSection data={data} />
      <LabsSection labs={data.labs} index={data.index} />
      <div className="grid grid-cols-1 gap-5 lg:grid-cols-3">
        <SectionCard title="Imaging observations">
          <ObservationGroupList groups={data.imaging} index={data.index} />
        </SectionCard>
        <SectionCard title="Assessments">
          <ObservationGroupList groups={data.assessments} index={data.index} />
        </SectionCard>
        <SectionCard title="Vitals & social history">
          <ObservationGroupList groups={data.vitals} index={data.index} />
        </SectionCard>
        <SectionCard title="Other conditions">
          <OtherConditionsList conditions={data.otherConditions} absent={data.comorbiditiesAbsent} notListed={data.comorbiditiesNotListed} />
        </SectionCard>
        <SectionCard title="Allergies">
          <AllergyList allergies={data.allergies} />
        </SectionCard>
      </div>
      <ReportsSection groups={data.reportGroups} />
      <OtherResources groups={data.otherResources} />
    </div>
  );
}
