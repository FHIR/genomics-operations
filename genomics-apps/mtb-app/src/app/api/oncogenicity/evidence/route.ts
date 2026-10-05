import { NextRequest, NextResponse } from 'next/server';
import { fetchEvidenceSummary } from '@/lib/oncogenicity';

export async function POST(request: NextRequest) {
    try {
        const body = await request.json() as { variant?: string; tumorType?: string };

        if (!body.variant) {
            return NextResponse.json({ error: 'variant is required' }, { status: 400 });
        }

        const evidence = await fetchEvidenceSummary(body.variant, body.tumorType);

        return NextResponse.json({ evidence });
    } catch (error) {
        const message = error instanceof Error ? error.message : 'Unknown error';
        return NextResponse.json({ error: message }, { status: 500 });
    }
}
