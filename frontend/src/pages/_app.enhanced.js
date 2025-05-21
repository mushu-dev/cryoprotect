import '../styles/globals.css';
import { CircuitBreakerProvider } from '../components/circuit-breaker';
import NavigationHeader from '../components/navigation-header-simple.jsx';

/**
 * Enhanced App component with CircuitBreakerProvider
 * 
 * This is a sample implementation showing how to integrate the
 * circuit breaker components into your app.
 */
function MyApp({ Component, pageProps }) {
  return (
    <CircuitBreakerProvider
      initialCircuits={['api', 'database', 'default']}
      refreshInterval={1000}
    >
      <div className="flex min-h-screen flex-col">
        <NavigationHeader />
        <main className="flex-1">
          <Component {...pageProps} />
        </main>
        <footer className="py-6 border-t">
          <div className="container mx-auto px-4">
            <div className="flex items-center justify-between">
              <div className="text-sm text-gray-500">
                &copy; {new Date().getFullYear()} CryoProtect
              </div>
            </div>
          </div>
        </footer>
      </div>
    </CircuitBreakerProvider>
  );
}

export default MyApp;