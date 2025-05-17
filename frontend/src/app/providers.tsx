'use client'

import { QueryClient, QueryClientProvider } from '@tanstack/react-query'
import { ReactNode, useState } from 'react'
import { ThemeProvider } from '@/components/theme-provider'
import { SessionProvider } from 'next-auth/react'
import { AnalyticsProvider } from '@/components/analytics/AnalyticsProvider'
import { AnalyticsConsent } from '@/components/analytics/AnalyticsConsent'

export function Providers({ children }: { children: ReactNode }) {
  const [queryClient] = useState(() => new QueryClient({
    defaultOptions: {
      queries: {
        staleTime: 60 * 1000, // 1 minute
        retry: 1,
        refetchOnWindowFocus: false,
      },
    },
  }))

  return (
    <ThemeProvider>
      <SessionProvider>
        <QueryClientProvider client={queryClient}>
          <AnalyticsProvider>
            {children}
            <AnalyticsConsent />
          </AnalyticsProvider>
        </QueryClientProvider>
      </SessionProvider>
    </ThemeProvider>
  )
}