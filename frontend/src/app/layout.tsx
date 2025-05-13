import type { Metadata } from 'next'
import { Inter } from 'next/font/google'
import '@/styles/globals.css'
import { Providers } from './providers'
import { NavigationHeader } from '@/components/navigation-header'
import { Footer } from '@/components/footer'

const inter = Inter({ 
  subsets: ['latin'],
  variable: '--font-sans',
})

export const metadata: Metadata = {
  title: 'CryoProtect - Cryoprotectant Analysis Platform',
  description: 'Analyze and optimize cryoprotectant molecules and mixtures',
  keywords: 'cryoprotectant, molecular analysis, scientific research, mixtures, visualization',
  authors: [{ name: 'CryoProtect Team' }],
}

export default function RootLayout({
  children,
}: {
  children: React.ReactNode
}) {
  return (
    <html lang="en" suppressHydrationWarning>
      <body className={`${inter.variable} font-sans min-h-screen flex flex-col`}>
        <Providers>
          <NavigationHeader />
          <main className="flex-grow container mx-auto py-6 px-4">
            {children}
          </main>
          <Footer />
        </Providers>
      </body>
    </html>
  )
}