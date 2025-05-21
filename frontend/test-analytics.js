/**
 * Analytics validation script
 * This script tests the functionality of the updated analytics integration
 * 
 * Run with: node test-analytics.js
 */

console.log('🧪 Testing Netlify Analytics Implementation');

// Simulate browser environment
global.window = {
  location: {
    pathname: '/test-page',
    hostname: 'cryoprotect.netlify.app'
  },
  localStorage: {
    _data: {},
    getItem(key) {
      return this._data[key];
    },
    setItem(key, value) {
      this._data[key] = value;
    },
    removeItem(key) {
      delete this._data[key];
    }
  },
  netlifyIdentity: {
    on: (event, callback) => {
      if (event === 'init') {
        callback({ id: 'test-user-id' });
      }
    }
  },
  _nenl: {
    push: (args) => {
      console.info(`🔍 Netlify Analytics Event:`, args);
    }
  },
  addEventListener: () => {},
  removeEventListener: () => {},
  document: {
    createElement: () => ({
      async: true,
      src: '',
      addEventListener: () => {},
    }),
    head: {
      appendChild: () => {},
      removeChild: () => {},
    }
  }
};

global.document = {
  createElement: () => ({
    async: true,
    src: '',
    addEventListener: () => {},
  }),
  head: {
    appendChild: () => {},
    removeChild: () => {},
  },
  dispatchEvent: () => {}
};

global.console.originalLog = console.log;
global.console.log = (message, ...args) => {
  console.info(`🔍 LOG: ${message}`, ...args);
};

// Import our analytics modules
const path = require('path');
const fs = require('fs');

// Helper to evaluate JS files
function evaluateModule(filePath) {
  const content = fs.readFileSync(filePath, 'utf8');
  const Module = module.constructor;
  const m = new Module();
  m._compile(`
    const React = { 
      useEffect: (fn) => fn(), 
      useState: (initial) => [initial, () => {}], 
      createContext: () => ({}),
      useContext: () => ({
        isEnabled: true,
        trackEvent: (name, props) => console.log('[Analytics] Event:', name, props),
        trackPageView: (path) => console.log('[Analytics] Page view:', path),
        trackFeatureUsage: (feature, props) => console.log('[Analytics] Feature used:', feature, props),
        setEnabled: (enabled) => console.log('[Analytics] Set enabled:', enabled)
      })
    };
    const usePathname = () => '/test-page';
    const useSearchParams = () => ({toString: () => ''});
    exports.useEffect = React.useEffect;
    exports.useState = React.useState;
    exports.createContext = React.createContext;
    exports.useContext = React.useContext;
    exports.usePathname = usePathname;
    exports.useSearchParams = useSearchParams;
    ${content}
  `, filePath);
  return m.exports;
}

// Set required environment variables
process.env.NEXT_PUBLIC_NETLIFY = 'true';
process.env.NODE_ENV = 'development';

// Test the netlify-analytics.js file
console.log('\n📝 Testing netlify-analytics.js');
try {
  // Mock necessary dependencies
  const analytics = require('./src/app/netlify-analytics');
  
  console.log('✅ Successfully loaded netlify-analytics.js');
  console.log('📊 Available exports:', Object.keys(analytics));
  
  // Test the trackEvent function
  if (typeof analytics.trackEvent === 'function') {
    analytics.trackEvent('test_event', { property: 'value' });
    console.log('✅ trackEvent executed successfully');
  } else {
    console.log('❌ trackEvent function not found');
  }
} catch (error) {
  console.error('❌ Error testing netlify-analytics.js:', error.message);
  console.error(error);
}

// Test the analytics hook
console.log('\n📝 Testing useAnalytics hook');
try {
  // Mock necessary dependencies
  const { useAnalytics } = require('./src/hooks/useAnalytics');
  
  // Create a mock hook environment
  const mockHook = () => {
    const result = {};
    // Mock the hook by calling it directly
    const hookResult = useAnalytics();
    // Copy properties
    Object.keys(hookResult).forEach(key => {
      result[key] = hookResult[key];
      // If it's a function, execute it to see if it works
      if (typeof hookResult[key] === 'function') {
        try {
          if (key === 'trackPageView') {
            hookResult[key]('/test-page');
          } else if (key === 'trackFeatureUsage') {
            hookResult[key]('test-feature', { detail: 'test' });
          } else if (key === 'trackError') {
            hookResult[key](new Error('test error'), 'test-source');
          } else if (key === 'trackSearch') {
            hookResult[key]('test query', 5);
          } else if (key === 'trackAction') {
            hookResult[key]('test action', 'test category', { detail: 'test' });
          }
        } catch (e) {
          console.log(`❌ Error executing ${key}:`, e.message);
        }
      }
    });
    return result;
  };
  
  const hookResult = mockHook();
  console.log('✅ useAnalytics hook executed successfully');
  console.log('📊 Hook exports:', Object.keys(hookResult));
  
} catch (error) {
  console.error('❌ Error testing useAnalytics hook:', error.message);
  console.error(error);
}

// Test the AnalyticsProvider component
console.log('\n📝 Testing AnalyticsProvider component');
try {
  // Mock React hooks more thoroughly for the provider
  const ReactMock = {
    createContext: (defaultValue) => ({
      Provider: ({ children, value }) => children,
      Consumer: ({ children }) => children(defaultValue),
      displayName: 'MockContext',
    }),
    useState: (initial) => {
      let state = initial;
      const setState = (newVal) => {
        state = typeof newVal === 'function' ? newVal(state) : newVal;
      };
      return [state, setState];
    },
    useContext: (ctx) => ({}),
    useEffect: (effect, deps) => effect(),
  };
  
  // Override module.exports before requiring
  const Module = module.constructor;
  Module.prototype._orig_compile = Module.prototype._compile;
  Module.prototype._compile = function(content, filename) {
    if (filename.includes('AnalyticsProvider.tsx')) {
      content = `
        const React = {
          createContext: (val) => ({
            Provider: ({ children }) => children,
            Consumer: ({ children }) => children(val)
          }),
          useState: (initial) => [initial, () => {}],
          useContext: () => ({}),
          useEffect: (fn) => fn()
        };
        const usePathname = () => '/test';
        const useSearchParams = () => ({ toString: () => '' });
        ${content.replace(/import[^;]*;/g, '').replace('export', 'module.exports =')}
      `;
    }
    return this._orig_compile(content, filename);
  };
  
  // Now let's test that we can parse the files, even if we can't fully execute them
  const AnalyticsProviderPath = path.join(__dirname, 'src/components/analytics/AnalyticsProvider.tsx');
  const AnalyticsConsentPath = path.join(__dirname, 'src/components/analytics/AnalyticsConsent.tsx');
  
  console.log('✅ Testing AnalyticsProvider and AnalyticsConsent parsing...');
  
  if (fs.existsSync(AnalyticsProviderPath)) {
    console.log('✅ AnalyticsProvider file exists');
  } else {
    console.log('❌ AnalyticsProvider file not found');
  }
  
  if (fs.existsSync(AnalyticsConsentPath)) {
    console.log('✅ AnalyticsConsent file exists');
  } else {
    console.log('❌ AnalyticsConsent file not found');
  }
  
  // Parse and compile just to test for syntax errors
  try {
    require('@/components/analytics/AnalyticsProvider');
    console.log('✅ AnalyticsProvider can be parsed');
  } catch (error) {
    console.log(`❌ Error parsing AnalyticsProvider: ${error.message.substring(0, 100)}...`);
  }
  
  try {
    require('@/components/analytics/AnalyticsConsent');
    console.log('✅ AnalyticsConsent can be parsed');
  } catch (error) {
    console.log(`❌ Error parsing AnalyticsConsent: ${error.message.substring(0, 100)}...`);
  }
  
  // Restore original compile function
  Module.prototype._compile = Module.prototype._orig_compile;
  delete Module.prototype._orig_compile;
  
} catch (error) {
  console.error('❌ Error testing AnalyticsProvider:', error.message);
}

// Test environment variable configuration
console.log('\n📝 Testing environment variable configuration');
try {
  // Check required environment variables
  const requiredVars = [
    'NEXT_PUBLIC_NETLIFY',
    'NEXT_PUBLIC_API_URL',
    'NEXT_PUBLIC_ENVIRONMENT'
  ];
  
  const netlifyToml = fs.readFileSync(path.join(__dirname, '..', 'netlify.toml'), 'utf8');
  const missingVars = [];
  
  requiredVars.forEach(varName => {
    if (!netlifyToml.includes(varName)) {
      missingVars.push(varName);
    }
  });
  
  if (missingVars.length > 0) {
    console.log(`⚠️ Missing environment variables in netlify.toml: ${missingVars.join(', ')}`);
  } else {
    console.log('✅ All required environment variables found in netlify.toml');
  }
  
  // Check .env.production as well
  const envProdPath = path.join(__dirname, '.env.production');
  if (fs.existsSync(envProdPath)) {
    const envProd = fs.readFileSync(envProdPath, 'utf8');
    const missingEnvVars = [];
    
    requiredVars.forEach(varName => {
      if (!envProd.includes(varName)) {
        missingEnvVars.push(varName);
      }
    });
    
    if (missingEnvVars.length > 0) {
      console.log(`⚠️ Missing environment variables in .env.production: ${missingEnvVars.join(', ')}`);
    } else {
      console.log('✅ All required environment variables found in .env.production');
    }
  } else {
    console.log('⚠️ .env.production file not found');
  }
  
} catch (error) {
  console.error('❌ Error testing environment variables:', error.message);
}

// Validate Content-Security-Policy
console.log('\n📝 Testing Content-Security-Policy');
try {
  const netlifyToml = fs.readFileSync(path.join(__dirname, '..', 'netlify.toml'), 'utf8');
  
  // Check for identity.netlify.com in CSP
  const netlifyIdentityDomain = 'identity.netlify.com';
  
  if (netlifyToml.includes(netlifyIdentityDomain)) {
    console.log(`✅ ${netlifyIdentityDomain} found in Content-Security-Policy`);
  } else {
    console.log(`⚠️ ${netlifyIdentityDomain} missing in Content-Security-Policy`);
  }
  
} catch (error) {
  console.error('❌ Error testing Content-Security-Policy:', error.message);
}

// Check for analytics-test page
console.log('\n📝 Testing analytics test page');
try {
  const analyticsTestPath = path.join(__dirname, 'src/pages/analytics-test.js');
  
  if (fs.existsSync(analyticsTestPath)) {
    console.log('✅ Analytics test page exists');
    
    // Check if the page imports the required components
    const content = fs.readFileSync(analyticsTestPath, 'utf8');
    
    if (content.includes('useAnalyticsContext')) {
      console.log('✅ Page imports analytics context');
    } else {
      console.log('⚠️ Page does not import analytics context');
    }
    
    if (content.includes('AnalyticsToggle')) {
      console.log('✅ Page includes analytics toggle component');
    } else {
      console.log('⚠️ Page does not include analytics toggle component');
    }
    
  } else {
    console.log('❌ Analytics test page not found');
  }
  
} catch (error) {
  console.error('❌ Error checking analytics test page:', error.message);
}

console.log('\n🏁 Analytics validation complete!');