'use client';

import React from 'react';
import * as Tooltip from '@radix-ui/react-tooltip';
import { Info } from 'lucide-react';

interface TooltipIconProps {
  content: string;
}

export default function TooltipIcon({ content }: TooltipIconProps) {
  return (
    <Tooltip.Provider>
      <Tooltip.Root>
        <Tooltip.Trigger asChild>
          <span
            className="inline-flex items-center justify-center text-gray-400 hover:text-gray-600 ml-1 cursor-default"
            aria-label="More information"
          >
            <Info size={16} />
          </span>
        </Tooltip.Trigger>
        <Tooltip.Portal>
          <Tooltip.Content
            side="top"
            sideOffset={4}
            className="bg-white border border-gray-300 rounded-lg shadow-md px-3 py-2 text-sm text-gray-700 max-w-sm z-50 animate-in fade-in-0 zoom-in-95 data-[state=closed]:animate-out data-[state=closed]:fade-out-0 data-[state=closed]:zoom-out-95"
          >
            {content}
            <Tooltip.Arrow className="fill-white stroke-gray-300" />
          </Tooltip.Content>
        </Tooltip.Portal>
      </Tooltip.Root>
    </Tooltip.Provider>
  );
}