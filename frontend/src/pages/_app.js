import '../styles/globals.css';
import NavigationHeader from '../components/navigation-header';

function MyApp({ Component, pageProps }) {
  return (
    <div className="flex min-h-screen flex-col">
      <NavigationHeader />
      <main className="flex-grow">
        <Component {...pageProps} />
      </main>
      <footer className="py-6 border-t">
        <div className="container mx-auto text-center text-muted-foreground">
          <p>© 2025 CryoProtect. All rights reserved.</p>
        </div>
      </footer>
    </div>
  );
}

export default MyApp;