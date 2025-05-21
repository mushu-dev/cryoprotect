import React from 'react';
import Head from 'next/head';
import Link from 'next/link';
import { 
  ArrowUpRight, 
  Beaker, 
  TrendingUp, 
  Thermometer, 
  Users, 
  Database,
  Flask,
  CircleCheck,
  AlertCircle
} from 'lucide-react';

import DashboardLayout from '../components/layouts/dashboard-layout';
import { 
  Card, 
  CardContent, 
  CardDescription, 
  CardFooter, 
  CardHeader, 
  CardTitle 
} from '../components/ui/card.jsx';
import { Separator } from '../components/ui/separator.jsx';

// Example stats
const stats = [
  {
    title: 'Total Molecules',
    value: '1,862',
    change: '+12%',
    trend: 'up',
    icon: <Database className="h-5 w-5" />,
    color: 'blue'
  },
  {
    title: 'Experiments',
    value: '284',
    change: '+18%',
    trend: 'up',
    icon: <Flask className="h-5 w-5" />,
    color: 'purple'
  },
  {
    title: 'Success Rate',
    value: '76%',
    change: '+5%',
    trend: 'up',
    icon: <CircleCheck className="h-5 w-5" />,
    color: 'green'
  },
  {
    title: 'Alerts',
    value: '3',
    change: '-2',
    trend: 'down',
    icon: <AlertCircle className="h-5 w-5" />,
    color: 'amber'
  }
];

// Example recent experiments
const recentExperiments = [
  {
    id: 'exp-001',
    title: 'High-Content Screening for DMSO Alternatives',
    status: 'Completed',
    date: '2025-05-15',
    success: true
  },
  {
    id: 'exp-002',
    title: 'Glycerol Efficacy in Stem Cell Cryopreservation',
    status: 'In Progress',
    date: '2025-05-18',
    success: null
  },
  {
    id: 'exp-003',
    title: 'Vitrification Protocol Optimization',
    status: 'Completed',
    date: '2025-05-14',
    success: true
  },
  {
    id: 'exp-004',
    title: 'Trehalose Concentration Series',
    status: 'Pending Analysis',
    date: '2025-05-13',
    success: null
  },
  {
    id: 'exp-005',
    title: 'Low-Toxicity Mixture Screening',
    status: 'Failed',
    date: '2025-05-10',
    success: false
  }
];

// Example top performing cryoprotectants
const topCryoprotectants = [
  {
    id: 'cryo-001',
    name: 'Trehalose',
    effectiveness: 94,
    usageCount: 78
  },
  {
    id: 'cryo-002',
    name: 'Dimethyl Sulfoxide',
    effectiveness: 92,
    usageCount: 156
  },
  {
    id: 'cryo-003',
    name: 'Glycerol',
    effectiveness: 89,
    usageCount: 112
  },
  {
    id: 'cryo-004',
    name: 'Propylene Glycol',
    effectiveness: 86,
    usageCount: 64
  },
  {
    id: 'cryo-005',
    name: 'Ethylene Glycol',
    effectiveness: 85,
    usageCount: 92
  }
];

export default function Dashboard() {
  return (
    <DashboardLayout>
      <Head>
        <title>Dashboard | CryoProtect</title>
        <meta name="description" content="CryoProtect dashboard with key metrics and performance indicators" />
      </Head>

      <div className="flex flex-col gap-8">
        <div>
          <h1 className="text-3xl font-bold tracking-tight">Dashboard</h1>
          <p className="text-muted-foreground mt-2">
            Welcome to CryoProtect. Here's an overview of your cryoprotectant experiments and data.
          </p>
        </div>

        {/* Stats Cards */}
        <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-4">
          {stats.map((stat, i) => (
            <Card key={i} className={`border-l-4 border-l-${stat.color}-500`}>
              <CardHeader className="flex flex-row items-center justify-between pb-2">
                <CardTitle className="text-sm font-medium">
                  {stat.title}
                </CardTitle>
                <div className={`p-2 rounded-full bg-${stat.color}-100`}>
                  {stat.icon}
                </div>
              </CardHeader>
              <CardContent>
                <div className="text-2xl font-bold">{stat.value}</div>
                <p className={`text-xs ${stat.trend === 'up' ? 'text-green-500' : 'text-red-500'} flex items-center`}>
                  {stat.change}
                  {stat.trend === 'up' ? (
                    <TrendingUp className="ml-1 h-3 w-3" />
                  ) : (
                    <ArrowUpRight className="ml-1 h-3 w-3 transform rotate-90" />
                  )}
                </p>
              </CardContent>
            </Card>
          ))}
        </div>

        <div className="grid gap-4 md:grid-cols-2">
          {/* Recent Experiments */}
          <Card className="col-span-1">
            <CardHeader>
              <CardTitle>Recent Experiments</CardTitle>
              <CardDescription>
                Your latest cryopreservation experiments
              </CardDescription>
            </CardHeader>
            <CardContent>
              <div className="space-y-4">
                {recentExperiments.map((experiment) => (
                  <div key={experiment.id} className="flex items-center justify-between">
                    <div>
                      <Link href={`/experiments/${experiment.id}`}>
                        <a className="font-medium hover:underline text-primary">
                          {experiment.title}
                        </a>
                      </Link>
                      <div className="text-sm text-muted-foreground">
                        {experiment.date}
                      </div>
                    </div>
                    <div className="flex items-center">
                      <span 
                        className={`inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-medium ${
                          experiment.status === 'Completed' ? 'bg-green-100 text-green-800' : 
                          experiment.status === 'Failed' ? 'bg-red-100 text-red-800' : 
                          'bg-blue-100 text-blue-800'
                        }`}
                      >
                        {experiment.status}
                      </span>
                    </div>
                  </div>
                ))}
              </div>
            </CardContent>
            <CardFooter>
              <Link href="/experiments">
                <a className="text-primary text-sm hover:underline">
                  View all experiments
                </a>
              </Link>
            </CardFooter>
          </Card>

          {/* Top Cryoprotectants */}
          <Card className="col-span-1">
            <CardHeader>
              <CardTitle>Top Performing Cryoprotectants</CardTitle>
              <CardDescription>
                Highest success rate based on experimental data
              </CardDescription>
            </CardHeader>
            <CardContent>
              <div className="space-y-4">
                {topCryoprotectants.map((cryo) => (
                  <div key={cryo.id} className="flex items-center justify-between">
                    <div className="flex items-center gap-2">
                      <Beaker className="h-4 w-4 text-primary" />
                      <Link href={`/molecules/${cryo.id}`}>
                        <a className="font-medium hover:underline text-primary">
                          {cryo.name}
                        </a>
                      </Link>
                    </div>
                    <div className="flex items-center gap-4">
                      <div className="text-sm">
                        <span className="text-muted-foreground">Used in: </span>
                        <span className="font-medium">{cryo.usageCount}</span>
                      </div>
                      <div>
                        <span 
                          className={`inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-medium ${
                            cryo.effectiveness > 90 ? 'bg-green-100 text-green-800' :
                            cryo.effectiveness > 80 ? 'bg-blue-100 text-blue-800' :
                            'bg-amber-100 text-amber-800'
                          }`}
                        >
                          {cryo.effectiveness}% Effective
                        </span>
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            </CardContent>
            <CardFooter>
              <Link href="/molecules">
                <a className="text-primary text-sm hover:underline">
                  View all molecules
                </a>
              </Link>
            </CardFooter>
          </Card>
        </div>

        {/* Quick Actions */}
        <Card>
          <CardHeader>
            <CardTitle>Quick Actions</CardTitle>
            <CardDescription>Common tasks and actions</CardDescription>
          </CardHeader>
          <CardContent>
            <div className="grid gap-4 grid-cols-1 sm:grid-cols-2 md:grid-cols-3">
              <Link href="/experiments/create">
                <a className="flex flex-col items-center justify-center p-4 bg-primary/5 hover:bg-primary/10 rounded-lg transition-colors">
                  <Flask className="h-8 w-8 text-primary mb-2" />
                  <span className="font-medium">New Experiment</span>
                </a>
              </Link>
              <Link href="/protocols/create">
                <a className="flex flex-col items-center justify-center p-4 bg-primary/5 hover:bg-primary/10 rounded-lg transition-colors">
                  <Thermometer className="h-8 w-8 text-primary mb-2" />
                  <span className="font-medium">Create Protocol</span>
                </a>
              </Link>
              <Link href="/molecules">
                <a className="flex flex-col items-center justify-center p-4 bg-primary/5 hover:bg-primary/10 rounded-lg transition-colors">
                  <Database className="h-8 w-8 text-primary mb-2" />
                  <span className="font-medium">Search Molecules</span>
                </a>
              </Link>
            </div>
          </CardContent>
        </Card>
      </div>
    </DashboardLayout>
  );
}