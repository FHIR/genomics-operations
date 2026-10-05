import { NextRequest, NextResponse } from 'next/server';
import { fetchPredictionObservation } from '@/lib/oncogenicity';

export async function POST(request: NextRequest) {
    try {
        const body = await request.json() as { variant?: string; tumorType?: string };

        if (!body.variant) {
            return NextResponse.json({ error: 'variant is required' }, { status: 400 });
        }

        const observation = await fetchPredictionObservation(body.variant, body.tumorType);

        return NextResponse.json({ observation });
    } catch (error) {
        const message = error instanceof Error ? error.message : 'Unknown error';
        return NextResponse.json({ error: message }, { status: 500 });
    }
}
