#!/bin/bash
set -e

echo "🔄 Simplifying App Router for build..."

# Backup App Router layout.tsx
if [ -f "src/app/layout.tsx" ]; then
  mkdir -p .app-router-backup
  cp src/app/layout.tsx .app-router-backup/
  
  # Create a simplified version
  cat > src/app/layout.tsx << EOL
import React from 'react';
import '../styles/globals.css';

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en">
      <body>
        <div className="flex min-h-screen flex-col bg-background">
          <header className="border-b py-4">
            <div className="container mx-auto">
              <h1 className="text-2xl font-bold">CryoProtect</h1>
            </div>
          </header>
          <main className="flex-grow">
            {children}
          </main>
          <footer className="border-t py-4">
            <div className="container mx-auto">
              <p className="text-sm text-gray-500">© 2025 CryoProtect</p>
            </div>
          </footer>
        </div>
      </body>
    </html>
  );
}
EOL
fi

echo "✅ App Router simplified for build"