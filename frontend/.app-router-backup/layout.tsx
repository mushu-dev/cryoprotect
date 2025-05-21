import React from 'react';
import '../styles/globals.css';
import { ThemeProvider } from '../components/theme-provider';
import NavigationHeader from '../components/navigation-header-updated';
import { Footer } from '../components/footer';
import { ConvexClientProvider } from '../convex/ConvexClientProvider';
import { Metadata, Viewport } from 'next';

export const metadata: Metadata = {
  title: {
    default: 'CryoProtect',
    template: '%s | CryoProtect',
  },
  description: 'A comprehensive platform for cryoprotectant analysis and experiment management',
  keywords: ['cryopreservation', 'cryoprotectant', 'science', 'molecules', 'experiments', 'protocols'],
  authors: [{ name: 'CryoProtect Team' }],
  creator: 'CryoProtect',
  openGraph: {
    type: 'website',
    locale: 'en_US',
    url: 'https://cryoprotect.app',
    title: 'CryoProtect',
    description: 'A comprehensive platform for cryoprotectant analysis and experiment management',
    siteName: 'CryoProtect',
  },
  twitter: {
    card: 'summary_large_image',
    title: 'CryoProtect',
    description: 'A comprehensive platform for cryoprotectant analysis and experiment management',
  },
  icons: {
    icon: '/favicon.ico',
    shortcut: '/favicon-16x16.png',
    apple: '/apple-touch-icon.png',
  },
  manifest: '/site.webmanifest',
};

export const viewport: Viewport = {
  themeColor: [
    { media: '(prefers-color-scheme: light)', color: '#ffffff' },
    { media: '(prefers-color-scheme: dark)', color: '#1e293b' },
  ],
  width: 'device-width',
  initialScale: 1,
  maximumScale: 5,
  userScalable: true,
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en" suppressHydrationWarning>
      <body className="antialiased">
        <ConvexClientProvider>
          <ThemeProvider attribute="class" defaultTheme="system" enableSystem>
            <div className="flex min-h-screen flex-col bg-background">
              <NavigationHeader />
              <main className="flex-grow">
                {children}
              </main>
              <Footer />
            </div>
          </ThemeProvider>
        </ConvexClientProvider>
      </body>
    </html>
  );
}