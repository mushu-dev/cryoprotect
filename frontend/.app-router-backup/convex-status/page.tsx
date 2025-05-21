import { ConvexStatus } from '../../components/convex-status'

export default function ConvexStatusPage() {
  return (
    <div className="container mx-auto py-12">
      <h1 className="text-3xl font-bold mb-8 text-center">Convex Database Status</h1>
      <p className="text-center mb-6 text-gray-600 dark:text-gray-400">
        This page allows you to check if your Convex database connection is working properly.
      </p>
      
      <div className="max-w-md mx-auto mt-8">
        <ConvexStatus />
      </div>
      
      <div className="mt-12 bg-gray-50 dark:bg-gray-800 rounded-lg p-6 max-w-2xl mx-auto">
        <h2 className="text-xl font-semibold mb-4">About Convex Integration</h2>
        <p className="mb-4">
          Convex is now the primary database for CryoProtect, offering:
        </p>
        <ul className="list-disc pl-6 space-y-2 mb-4">
          <li>Real-time data synchronization</li>
          <li>Automatic caching and query optimization</li>
          <li>TypeScript integration for type safety</li>
          <li>Scalable database infrastructure</li>
          <li>Simplified authentication flow</li>
        </ul>
        <p>
          Convex is configured via environment variables and properly integrated with the application.
          If you see "Connected" status above, your Convex database is working correctly.
        </p>
      </div>
    </div>
  )
}