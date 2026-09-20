import { NextRequest, NextResponse } from 'next/server';
import { convertSpdiToHgvs, fetchEvidenceSummary } from '@/lib/oncogenicity';

export async function POST(request: NextRequest) {
    try {
        const body = await request.json() as { spdi?: string };

        if (!body.spdi) {
            return NextResponse.json({ error: 'spdi is required' }, { status: 400 });
        }

        const hgvs = await convertSpdiToHgvs(body.spdi);
        const evidence = await fetchEvidenceSummary(hgvs);

        return NextResponse.json({ hgvs, evidence });
    } catch (error) {
        const message = error instanceof Error ? error.message : 'Unknown error';
        return NextResponse.json({ error: message }, { status: 500 });
    }
}
