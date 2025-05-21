/**
 * Home Page
 * Main entry point for the application with links to molecules and mixtures
 */
import React from 'react';
import Link from 'next/link';
import { Button } from '../components/ui/button';
import { Card, CardContent, CardHeader, CardFooter, CardTitle } from '../components/ui/card';

export default function HomePage() {
  return (
    <div className="container mx-auto py-16">
      <div className="text-center mb-12">
        <h1 className="text-4xl font-bold mb-4">CryoProtect</h1>
        <p className="text-xl text-muted-foreground max-w-2xl mx-auto">
          A comprehensive platform for exploring, analyzing, and designing cryoprotectants and their mixtures.
        </p>
      </div>
      
      <div className="grid grid-cols-1 md:grid-cols-2 gap-8 max-w-4xl mx-auto">
        {/* Molecules Card */}
        <Card className="flex flex-col h-full">
          <CardHeader>
            <CardTitle>Molecules</CardTitle>
          </CardHeader>
          <CardContent className="flex-grow">
            <p>
              Browse our database of molecules, including cryoprotectants and related compounds. 
              View molecular structures, properties, and detailed information.
            </p>
          </CardContent>
          <CardFooter>
            <Button asChild className="w-full">
              <Link href="/molecules">
                View Molecules
              </Link>
            </Button>
          </CardFooter>
        </Card>
        
        {/* Mixtures Card */}
        <Card className="flex flex-col h-full">
          <CardHeader>
            <CardTitle>Mixtures</CardTitle>
          </CardHeader>
          <CardContent className="flex-grow">
            <p>
              Explore cryoprotectant mixtures and formulations. View composition charts, 
              cryoprotection scores, and analyze mixture properties.
            </p>
          </CardContent>
          <CardFooter>
            <Button asChild className="w-full">
              <Link href="/mixtures">
                View Mixtures
              </Link>
            </Button>
          </CardFooter>
        </Card>
      </div>
      
      <div className="mt-16 text-center text-sm text-muted-foreground">
        <p>
          Built with Convex + Next.js
        </p>
        <p className="mt-2">
          © 2025 CryoProtect. All rights reserved.
        </p>
      </div>
    </div>
  );
}