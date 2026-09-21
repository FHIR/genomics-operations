import { NextRequest, NextResponse } from 'next/server';
import { convertSpdiToHgvs, fetchPredictionObservation } from '@/lib/oncogenicity';

export async function POST(request: NextRequest) {
    try {
        const body = await request.json() as { spdi?: string; tumorType?: string };

        if (!body.spdi) {
            return NextResponse.json({ error: 'spdi is required' }, { status: 400 });
        }

        const hgvs = await convertSpdiToHgvs(body.spdi);
        const observation = await fetchPredictionObservation(hgvs, body.tumorType);

        return NextResponse.json({ hgvs, observation });
    } catch (error) {
        const message = error instanceof Error ? error.message : 'Unknown error';
        return NextResponse.json({ error: message }, { status: 500 });
    }
}
