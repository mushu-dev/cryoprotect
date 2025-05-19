import React from 'react';
import '../styles/globals.css';
import { ThemeProvider } from '../components/theme-provider';
import NavigationHeader from '../components/navigation-header';
import Footer from '../components/footer';

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en">
      <head>
        <title>CryoProtect</title>
        <meta name="description" content="A platform for cryoprotectant analysis" />
      </head>
      <body>
        <ThemeProvider attribute="class" defaultTheme="system" enableSystem>
          <div className="flex min-h-screen flex-col">
            <NavigationHeader />
            <main className="flex-grow">
              {children}
            </main>
            <Footer />
          </div>
        </ThemeProvider>
      </body>
    </html>
  );
}