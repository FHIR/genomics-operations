interface LoadingCellProps {
  isLoading?: boolean;
  error?: string;
  loadingMessage?: string;
  children: React.ReactNode;
}

export default function LoadingCell({
  isLoading,
  error,
  children,
  loadingMessage = "Loading..."
}: LoadingCellProps) {
  if (error) {
    return (
      <div className="p-3 text-red-400 text-sm">
        <div className="flex items-center gap-2">
          <svg xmlns="http://www.w3.org/2000/svg" className="h-4 w-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-2.5L13.732 4c-.77-.833-1.964-.833-2.732 0L4.082 16.5c-.77.833.192 2.5 1.732 2.5z" />
          </svg>
          <span>Error: {error}</span>
        </div>
      </div>
    );
  }

  if (isLoading) {
    return (
      <div className="p-3 text-gray-400 text-sm">
        <div className="flex items-center gap-2">
          <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-blue-500"></div>
          <span>{loadingMessage}</span>
        </div>
      </div>
    );
  }

  return <>{children}</>;
}