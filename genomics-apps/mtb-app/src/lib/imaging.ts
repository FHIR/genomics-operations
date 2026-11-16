import { DiagnosticReport, Endpoint, FhirResource, ImagingStudy } from '@/types/fhir';
import { conceptLabel } from './clinicalSummary';
import { reportDate, resolveReference, resourceKey } from './caseSummary';

// The NCI Imaging Data Commons (IDC) hosts the images; its viewers are embedded in the Imaging tab
const IDC_HOST = 'imaging.datacommons.cancer.gov';
const IDC_VIEWER_BASE = 'https://viewer.imaging.datacommons.cancer.gov';

export interface ImageEntry {
  key: string;
  study: ImagingStudy;
  studyUid?: string;
  // DICOM modality codes, e.g. ['CT', 'PT']
  modalities: string[];
  // Short type for the table: CT, PET/CT, MR, NM, Slide
  typeLabel: string;
  isSlide: boolean;
  slideCount: number;
  date?: string;
  description: string;
  // Reports that link this study, newest first
  reports: DiagnosticReport[];
  // The IDC DICOMweb service holding the images, when the study is in IDC
  dicomWebBase?: string;
  // The study date in the images themselves, when the record's note gives it (synthetic records shift their dates)
  imageDate?: string;
}

function studyUid(study: ImagingStudy) {
  const value = study.identifier?.find((identifier) => identifier.system === 'urn:dicom:uid')?.value;
  return value?.replace(/^urn:oid:/, '');
}

function modalityCodes(study: ImagingStudy) {
  const codes = [...(study.modality ?? []), ...(study.series ?? []).map((series) => series.modality).filter((coding) => coding !== undefined)]
    .map((coding) => coding.code)
    .filter((code): code is string => Boolean(code));
  return [...new Set(codes)];
}

function typeLabel(modalities: string[]) {
  if (modalities.includes('SM')) {
    return 'Slide';
  }
  if (modalities.includes('PT') && modalities.includes('CT')) {
    return 'PET/CT';
  }
  return modalities.join('/') || 'Image';
}

// Synthetic records state the real image date in the study's note, e.g. "The study date in the DICOM headers (shown by the IDC viewer) is 2001-12-29."
function noteImageDate(study: ImagingStudy) {
  for (const note of study.note ?? []) {
    const match = note.text?.match(/DICOM headers[^.]*? is (\d{4}-\d{2}-\d{2})/);
    if (match) {
      return match[1];
    }
  }
  return undefined;
}

function idcDicomWebBase(study: ImagingStudy, index: Map<string, FhirResource>) {
  const references = [...(study.endpoint ?? []), ...(study.series ?? []).flatMap((series) => series.endpoint ?? [])];
  for (const reference of references) {
    const endpoint = resolveReference(reference, index) as Endpoint | undefined;
    if (endpoint?.address?.includes(IDC_HOST)) {
      return endpoint.address.replace(/\/$/, '');
    }
  }
  return undefined;
}

/**
 * Every ImagingStudy in the record, newest first, with the reports that link it.
 */
export function buildImageCatalog(resources: FhirResource[]): ImageEntry[] {
  const index = new Map(resources.map((resource) => [resourceKey(resource), resource]));
  const reports = resources.filter((resource): resource is DiagnosticReport => resource.resourceType === 'DiagnosticReport'
    && resource.status !== 'entered-in-error');
  const reportsByStudy = new Map<string, DiagnosticReport[]>();
  for (const report of reports) {
    for (const reference of report.imagingStudy ?? []) {
      const study = resolveReference(reference, index);
      if (study) {
        const key = resourceKey(study);
        reportsByStudy.set(key, [...(reportsByStudy.get(key) ?? []), report]);
      }
    }
  }

  return resources
    .filter((resource): resource is ImagingStudy => resource.resourceType === 'ImagingStudy' && resource.status !== 'entered-in-error')
    .map((study) => {
      const key = resourceKey(study);
      const modalities = modalityCodes(study);
      const isSlide = modalities.includes('SM');
      return {
        key,
        study,
        studyUid: studyUid(study),
        modalities,
        typeLabel: typeLabel(modalities),
        isSlide,
        slideCount: isSlide ? (study.series ?? []).filter((series) => series.modality?.code === 'SM').length : 0,
        date: study.started?.slice(0, 10),
        description: study.description ?? (isSlide ? 'Pathology slides' : 'Imaging study'),
        reports: (reportsByStudy.get(key) ?? []).sort((left, right) => (reportDate(right) ?? '').localeCompare(reportDate(left) ?? '')),
        dicomWebBase: idcDicomWebBase(study, index),
        imageDate: noteImageDate(study),
      };
    })
    .sort((left, right) => (right.date ?? '').localeCompare(left.date ?? '') || left.description.localeCompare(right.description));
}

export function reportTitle(report: DiagnosticReport) {
  return conceptLabel(report.code) ?? 'Report';
}

/**
 * The image link's label for a report: "3 slides" when it links only slide studies, else "2 image studies".
 */
export function reportImagesLabel(entries: ImageEntry[]) {
  if (entries.length > 0 && entries.every((entry) => entry.isSlide)) {
    const slides = entries.reduce((total, entry) => total + Math.max(entry.slideCount, 1), 0);
    return `${slides} ${slides === 1 ? 'slide' : 'slides'}`;
  }
  return `${entries.length} ${entries.length === 1 ? 'image study' : 'image studies'}`;
}

/**
 * Each report's linked images, keyed by the report's resource key.
 */
export function imagesByReport(catalog: ImageEntry[]) {
  const byReport = new Map<string, ImageEntry[]>();
  for (const entry of catalog) {
    for (const report of entry.reports) {
      const key = resourceKey(report);
      byReport.set(key, [...(byReport.get(key) ?? []), entry]);
    }
  }
  return byReport;
}

/**
 * The IDC viewer URL for a study: the slide viewer (SLIM) for slides, opening the given slide when known, else
 * the radiology viewer (OHIF) on the study's first series. Undefined when the study isn't in IDC.
 */
export function idcViewerUrl(entry: ImageEntry, slideSeriesUid?: string) {
  if (!entry.dicomWebBase || !entry.studyUid) {
    return undefined;
  }
  if (entry.isSlide) {
    return `${IDC_VIEWER_BASE}/slim/studies/${entry.studyUid}${slideSeriesUid ? `/series/${slideSeriesUid}` : ''}`;
  }
  const firstSeries = entry.study.series?.[0]?.uid;
  return `${IDC_VIEWER_BASE}/v3/viewer/?StudyInstanceUIDs=${entry.studyUid}${firstSeries ? `&initialSeriesInstanceUID=${firstSeries}` : ''}`;
}
