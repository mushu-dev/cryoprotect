import React, { useState } from 'react';
import Head from 'next/head';
import Link from 'next/link';

export default function Layout({ children, title = 'CryoProtect' }) {
  const [isMobileMenuOpen, setIsMobileMenuOpen] = useState(false);

  const toggleMobileMenu = () => {
    setIsMobileMenuOpen(!isMobileMenuOpen);
    
    // When opening the menu, scroll to top to ensure it's visible
    if (!isMobileMenuOpen) {
      window.scrollTo({
        top: 0,
        behavior: 'smooth'
      });
    }
  };

  return (
    <>
      <Head>
        <title>{title} | CryoProtect</title>
        <meta name="viewport" content="width=device-width, initial-scale=1.0" />
        <meta name="description" content="Cryoprotectant research and analysis platform" />
        <link rel="icon" href="/favicon.ico" />
      </Head>

      <header className="header">
        <div className="header-container">
          <div className="header-top-row">
            <Link href="/">
              <span className="header-logo">CryoProtect</span>
            </Link>
            
            <button 
              className={`mobile-menu-toggle ${isMobileMenuOpen ? 'open' : ''}`} 
              onClick={toggleMobileMenu}
              aria-label="Toggle navigation menu"
            >
              <span className="toggle-bar"></span>
              <span className="toggle-bar"></span>
              <span className="toggle-bar"></span>
            </button>
            
            <div className="secondary-nav desktop-only">
              <Link href="/search" className="nav-link-search" title="Search">
                <svg width="20" height="20" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z"></path>
                </svg>
              </Link>
              <Link href="/convex" className="nav-link-secondary">
                Convex Demo
              </Link>
              <Link href="/connections" className="nav-link-secondary">
                Connections
              </Link>
              <Link href="/about" className="nav-link-secondary">
                About
              </Link>
            </div>
          </div>
          
          <nav className={`main-nav ${isMobileMenuOpen ? 'open' : ''}`}>
            <ul className="nav-links">
              <li>
                <Link href="/dashboard" className="nav-link" onClick={() => setIsMobileMenuOpen(false)}>
                  <span className="nav-icon">📊</span>
                  <span>Dashboard</span>
                </Link>
              </li>
              <li>
                <Link href="/molecules" className="nav-link" onClick={() => setIsMobileMenuOpen(false)}>
                  <span className="nav-icon">🧪</span>
                  <span>Molecules</span>
                </Link>
              </li>
              <li>
                <Link href="/mixtures" className="nav-link" onClick={() => setIsMobileMenuOpen(false)}>
                  <span className="nav-icon">🧬</span>
                  <span>Mixtures</span>
                </Link>
              </li>
              <li>
                <Link href="/protocols" className="nav-link" onClick={() => setIsMobileMenuOpen(false)}>
                  <span className="nav-icon">📋</span>
                  <span>Protocols</span>
                </Link>
              </li>
              <li>
                <Link href="/experiments" className="nav-link" onClick={() => setIsMobileMenuOpen(false)}>
                  <span className="nav-icon">🔬</span>
                  <span>Experiments</span>
                </Link>
              </li>
            </ul>
          </nav>
          
          <div className="secondary-nav mobile-only">
            <Link href="/search" className="nav-link-search" title="Search" onClick={() => setIsMobileMenuOpen(false)}>
              <svg width="20" height="20" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z"></path>
              </svg>
              <span>Search</span>
            </Link>
            <Link href="/convex" className="nav-link-secondary" onClick={() => setIsMobileMenuOpen(false)}>
              Convex Demo
            </Link>
            <Link href="/connections" className="nav-link-secondary" onClick={() => setIsMobileMenuOpen(false)}>
              Connections
            </Link>
            <Link href="/about" className="nav-link-secondary" onClick={() => setIsMobileMenuOpen(false)}>
              About
            </Link>
          </div>
        </div>
        
        <style jsx>{`
          .header-top-row {
            display: flex;
            justify-content: space-between;
            align-items: center;
            width: 100%;
          }
          
          .mobile-menu-toggle {
            display: none;
            flex-direction: column;
            justify-content: space-between;
            width: 30px;
            height: 21px;
            background: transparent;
            border: none;
            cursor: pointer;
            padding: 0;
            z-index: 10;
          }
          
          .toggle-bar {
            width: 100%;
            height: 3px;
            background-color: #333;
            border-radius: 2px;
            transition: all 0.3s;
          }
          
          .mobile-menu-toggle.open .toggle-bar:nth-child(1) {
            transform: translateY(9px) rotate(45deg);
          }
          
          .mobile-menu-toggle.open .toggle-bar:nth-child(2) {
            opacity: 0;
          }
          
          .mobile-menu-toggle.open .toggle-bar:nth-child(3) {
            transform: translateY(-9px) rotate(-45deg);
          }
          
          .desktop-only {
            display: flex;
          }
          
          .mobile-only {
            display: none;
          }
          
          @media (max-width: 768px) {
            .mobile-menu-toggle {
              display: flex;
            }
            
            .main-nav {
              overflow: hidden;
              max-height: 0;
              transition: max-height 0.3s ease-in-out;
              width: 100%;
            }
            
            .main-nav.open {
              max-height: 1000px;
            }
            
            .desktop-only {
              display: none;
            }
            
            .mobile-only {
              display: flex;
              flex-direction: column;
              width: 100%;
              margin-top: 10px;
            }
            
            .mobile-only .nav-link-search {
              width: 100%;
              height: auto;
              border-radius: 4px;
              padding: 8px 12px;
              justify-content: flex-start;
              gap: 8px;
            }
            
            .mobile-only .nav-link-secondary {
              padding: 8px 12px;
              margin: 4px 0;
            }
          }
        `}</style>
      </header>

      <main className="main">
        {children}
      </main>

      <footer className="footer">
        <div className="footer-content">
          <div className="footer-copyright">
            © {new Date().getFullYear()} CryoProtect - All rights reserved
          </div>
          <div className="footer-links">
            <Link href="/about">About</Link>
            <span className="footer-separator">|</span>
            <Link href="/">Home</Link>
          </div>
        </div>
      </footer>
    </>
  );
}