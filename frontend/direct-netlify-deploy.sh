#!/bin/bash
set -e

echo "🚀 Running direct Netlify deploy script..."

# Create a directory for static files
mkdir -p netlify-static-site/out

# Create HTML files
echo "📝 Creating static HTML files..."

# Create index.html
cat > netlify-static-site/out/index.html << EOL
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>CryoProtect - Cryoprotectant Analysis Platform</title>
  <meta name="description" content="Analyze and optimize cryoprotectant molecules and mixtures">
  <style>
    body {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans', 'Helvetica Neue', sans-serif;
      line-height: 1.6;
      color: #333;
      max-width: 800px;
      margin: 0 auto;
      padding: 2rem;
    }
    h1 {
      font-size: 2.5rem;
      margin-bottom: 1rem;
      color: #0070f3;
    }
    p {
      margin: 1rem 0;
    }
    .button {
      display: inline-block;
      background-color: #0070f3;
      color: white;
      padding: 0.5rem 1rem;
      text-decoration: none;
      border-radius: 0.25rem;
      font-weight: 500;
      margin-right: 0.5rem;
      margin-bottom: 0.5rem;
    }
    .button:hover {
      background-color: #0060df;
    }
    .container {
      margin: 2rem 0;
    }
    .card {
      border: 1px solid #eaeaea;
      border-radius: 0.5rem;
      padding: 1.5rem;
      margin: 1rem 0;
    }
    .card h2 {
      margin-top: 0;
      color: #0070f3;
    }
    .footer {
      margin-top: 3rem;
      padding-top: 1rem;
      border-top: 1px solid #eaeaea;
      font-size: 0.875rem;
      color: #666;
    }
    @media (max-width: 600px) {
      body {
        padding: 1rem;
      }
    }
  </style>
</head>
<body>
  <h1>CryoProtect</h1>
  <p>Welcome to the CryoProtect Cryoprotectant Analysis Platform</p>
  
  <div class="container">
    <a href="/molecules" class="button">Explore Molecules</a>
    <a href="/mixtures" class="button">View Mixtures</a>
    <a href="/properties" class="button">Search Properties</a>
  </div>
  
  <div class="card">
    <h2>Getting Started</h2>
    <p>
      CryoProtect helps researchers analyze and optimize cryoprotectant molecules and mixtures.
      Browse our extensive database of cryoprotectant molecules, create custom mixtures,
      and analyze their properties.
    </p>
  </div>
  
  <div class="card">
    <h2>Features</h2>
    <ul>
      <li>Comprehensive molecule database</li>
      <li>Mixture composition analysis</li>
      <li>Property predictions</li>
      <li>3D molecular visualization</li>
      <li>Data export capabilities</li>
    </ul>
  </div>
  
  <div class="footer">
    <p>
      CryoProtect - Advanced Cryoprotectant Analysis Platform
    </p>
  </div>
  
  <script>
    // Client-side navigation handler
    document.addEventListener('DOMContentLoaded', function() {
      // Add event listeners to navigation links
      document.querySelectorAll('a').forEach(function(link) {
        link.addEventListener('click', function(e) {
          // Get the href attribute
          var href = this.getAttribute('href');
          
          // Check if it's a relative link and not an external link
          if (href && href.indexOf('/') === 0 && href.indexOf('//') !== 0) {
            e.preventDefault();
            
            // Update browser history
            window.history.pushState({}, '', href);
            
            // You would normally load the new page content here via AJAX
            // For now, we'll just show a loading message
            document.body.innerHTML = '<div style="text-align: center; padding: 3rem;"><h1>Loading...</h1><p>Please wait while we load the page.</p><p><a href="/">Return Home</a></p></div>';
            
            // Redirect after a short delay to simulate loading
            setTimeout(function() {
              window.location.href = href;
            }, 500);
          }
        });
      });
      
      // Handle browser back/forward navigation
      window.addEventListener('popstate', function() {
        window.location.reload();
      });
    });
  </script>
</body>
</html>
EOL

# Create 404.html
cat > netlify-static-site/out/404.html << EOL
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>404 - Page Not Found | CryoProtect</title>
  <style>
    body {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans', 'Helvetica Neue', sans-serif;
      line-height: 1.6;
      color: #333;
      max-width: 600px;
      margin: 0 auto;
      padding: 2rem;
      text-align: center;
    }
    h1 {
      font-size: 2rem;
      margin-bottom: 1rem;
      color: #0070f3;
    }
    p {
      margin: 1rem 0;
    }
    .button {
      display: inline-block;
      background-color: #0070f3;
      color: white;
      padding: 0.5rem 1rem;
      text-decoration: none;
      border-radius: 0.25rem;
      font-weight: 500;
      margin-top: 1rem;
    }
    .button:hover {
      background-color: #0060df;
    }
  </style>
</head>
<body>
  <h1>404 - Page Not Found</h1>
  <p>The page you are looking for does not exist or has been moved.</p>
  <a href="/" class="button">Return to Home</a>
</body>
</html>
EOL

# Create molecules/[id].html
mkdir -p netlify-static-site/out/molecules
cat > netlify-static-site/out/molecules/index.html << EOL
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Molecules | CryoProtect</title>
  <style>
    body {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans', 'Helvetica Neue', sans-serif;
      line-height: 1.6;
      color: #333;
      max-width: 800px;
      margin: 0 auto;
      padding: 2rem;
    }
    h1 {
      font-size: 2rem;
      margin-bottom: 1rem;
      color: #0070f3;
    }
    .card {
      border: 1px solid #eaeaea;
      border-radius: 0.5rem;
      padding: 1.5rem;
      margin: 1rem 0;
    }
    .button {
      display: inline-block;
      background-color: #0070f3;
      color: white;
      padding: 0.5rem 1rem;
      text-decoration: none;
      border-radius: 0.25rem;
      font-weight: 500;
      margin-right: 0.5rem;
    }
    .home-link {
      display: inline-flex;
      align-items: center;
      color: #0070f3;
      text-decoration: none;
      margin-bottom: 1rem;
    }
    .molecule-list {
      list-style: none;
      padding: 0;
    }
    .molecule-item {
      padding: 0.75rem;
      border-bottom: 1px solid #eaeaea;
    }
    .molecule-link {
      color: #0070f3;
      text-decoration: none;
      display: block;
    }
    .molecule-link:hover {
      text-decoration: underline;
    }
  </style>
</head>
<body>
  <a href="/" class="home-link">
    <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M19 12H5M12 19l-7-7 7-7"/></svg>
    <span style="margin-left: 0.5rem;">Home</span>
  </a>
  
  <h1>Molecules</h1>
  <p>Browse our collection of cryoprotectant molecules.</p>
  
  <div class="card">
    <h2>Featured Molecules</h2>
    <ul class="molecule-list">
      <li class="molecule-item">
        <a href="/molecules/962" class="molecule-link">Glycerol (CID: 962)</a>
        <p>A common cryoprotectant used in various applications</p>
      </li>
      <li class="molecule-item">
        <a href="/molecules/176" class="molecule-link">DMSO - Dimethyl sulfoxide (CID: 176)</a>
        <p>Widely used in cell preservation and cryobiology</p>
      </li>
      <li class="molecule-item">
        <a href="/molecules/6276" class="molecule-link">Ethylene Glycol (CID: 6276)</a>
        <p>Effective cryoprotectant for various biological tissues</p>
      </li>
      <li class="molecule-item">
        <a href="/molecules/8857" class="molecule-link">Propylene Glycol (CID: 8857)</a>
        <p>Used in specific cryopreservation applications</p>
      </li>
    </ul>
    
    <div style="margin-top: 1.5rem;">
      <a href="/molecules/import" class="button">Import New Molecule</a>
    </div>
  </div>
  
  <p style="margin-top: 2rem; color: #666; font-size: 0.875rem;">
    For the full molecular database experience, please use the main application with JavaScript enabled.
  </p>

  <script>
    // Client-side navigation handler
    document.addEventListener('DOMContentLoaded', function() {
      // Add event listeners to navigation links
      document.querySelectorAll('a').forEach(function(link) {
        link.addEventListener('click', function(e) {
          // Get the href attribute
          var href = this.getAttribute('href');
          
          // Check if it's a relative link and not an external link
          if (href && href.indexOf('/') === 0 && href.indexOf('//') !== 0) {
            e.preventDefault();
            
            // Update browser history
            window.history.pushState({}, '', href);
            
            // You would normally load the new page content here via AJAX
            // For now, we'll just show a loading message
            document.body.innerHTML = '<div style="text-align: center; padding: 3rem;"><h1>Loading...</h1><p>Please wait while we load the page.</p><p><a href="/">Return Home</a></p></div>';
            
            // Redirect after a short delay to simulate loading
            setTimeout(function() {
              window.location.href = href;
            }, 500);
          }
        });
      });
      
      // Handle browser back/forward navigation
      window.addEventListener('popstate', function() {
        window.location.reload();
      });
    });
  </script>
</body>
</html>
EOL

# Create molecule detail page template
cat > netlify-static-site/out/molecules/detail.html << EOL
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Molecule Details | CryoProtect</title>
  <style>
    body {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans', 'Helvetica Neue', sans-serif;
      line-height: 1.6;
      color: #333;
      max-width: 800px;
      margin: 0 auto;
      padding: 2rem;
    }
    h1 {
      font-size: 2rem;
      margin-bottom: 0.5rem;
      color: #0070f3;
    }
    .home-link {
      display: inline-flex;
      align-items: center;
      color: #0070f3;
      text-decoration: none;
      margin-bottom: 1rem;
    }
    .card {
      border: 1px solid #eaeaea;
      border-radius: 0.5rem;
      padding: 1.5rem;
      margin: 1rem 0;
    }
    .badge {
      display: inline-block;
      background-color: rgba(0, 112, 243, 0.1);
      color: #0070f3;
      padding: 0.25rem 0.5rem;
      border-radius: 0.25rem;
      font-size: 0.875rem;
      border: 1px solid rgba(0, 112, 243, 0.2);
      margin-left: 0.5rem;
    }
    .dl-horizontal {
      display: grid;
      grid-template-columns: 1fr 2fr;
      row-gap: 0.75rem;
      column-gap: 1rem;
    }
    .dl-horizontal dt {
      font-weight: 500;
      color: #666;
    }
    .tabs {
      display: flex;
      border-bottom: 1px solid #eaeaea;
      margin-bottom: 1rem;
    }
    .tab {
      padding: 0.5rem 1rem;
      cursor: pointer;
      color: #666;
      border-bottom: 2px solid transparent;
    }
    .tab.active {
      color: #0070f3;
      border-bottom-color: #0070f3;
    }
    .tab-content {
      padding: 1rem 0;
    }
    .tab-content:not(.active) {
      display: none;
    }
    .loading {
      display: flex;
      flex-direction: column;
      align-items: center;
      justify-content: center;
      height: 300px;
      background-color: #f9f9f9;
      border-radius: 0.5rem;
    }
  </style>
</head>
<body>
  <a href="/molecules" class="home-link">
    <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M19 12H5M12 19l-7-7 7-7"/></svg>
    <span style="margin-left: 0.5rem;">Back to Molecules</span>
  </a>
  
  <h1>Molecule Details <span class="badge">Cryoprotectant</span></h1>
  <p>Loading molecule information...</p>
  
  <div class="card">
    <h2>Molecule Information</h2>
    <div class="loading">
      <svg width="38" height="38" viewBox="0 0 38 38" xmlns="http://www.w3.org/2000/svg" stroke="#0070f3">
        <g fill="none" fill-rule="evenodd">
          <g transform="translate(1 1)" stroke-width="2">
            <circle stroke-opacity=".5" cx="18" cy="18" r="18"/>
            <path d="M36 18c0-9.94-8.06-18-18-18">
              <animateTransform attributeName="transform" type="rotate" from="0 18 18" to="360 18 18" dur="1s" repeatCount="indefinite"/>
            </path>
          </g>
        </g>
      </svg>
      <p style="margin-top: 1rem;">Loading molecule data...</p>
    </div>
  </div>
  
  <div style="margin-top: 2rem;">
    <div class="tabs">
      <div class="tab active" data-tab="properties">Properties</div>
      <div class="tab" data-tab="visualization">Visualization</div>
      <div class="tab" data-tab="mixtures">Related Mixtures</div>
    </div>
    
    <div class="tab-content active" id="properties-tab">
      <div class="card">
        <div class="loading">
          <p>Loading property data...</p>
        </div>
      </div>
    </div>
    
    <div class="tab-content" id="visualization-tab">
      <div class="card">
        <div class="loading">
          <p>3D visualization requires JavaScript to be enabled</p>
        </div>
      </div>
    </div>
    
    <div class="tab-content" id="mixtures-tab">
      <div class="card">
        <div class="loading">
          <p>Loading related mixtures...</p>
        </div>
      </div>
    </div>
  </div>
  
  <p style="margin-top: 2rem; color: #666; font-size: 0.875rem;">
    For the full molecular visualization experience, please use the main application with JavaScript enabled.
  </p>

  <script>
    // Tab navigation
    document.querySelectorAll('.tab').forEach(function(tab) {
      tab.addEventListener('click', function() {
        // Remove active class from all tabs
        document.querySelectorAll('.tab').forEach(function(t) {
          t.classList.remove('active');
        });
        
        // Add active class to clicked tab
        this.classList.add('active');
        
        // Hide all tab contents
        document.querySelectorAll('.tab-content').forEach(function(content) {
          content.classList.remove('active');
        });
        
        // Show the selected tab content
        var tabId = this.getAttribute('data-tab');
        document.getElementById(tabId + '-tab').classList.add('active');
      });
    });
    
    // Client-side navigation handler
    document.addEventListener('DOMContentLoaded', function() {
      document.querySelectorAll('a').forEach(function(link) {
        link.addEventListener('click', function(e) {
          var href = this.getAttribute('href');
          
          if (href && href.indexOf('/') === 0 && href.indexOf('//') !== 0) {
            e.preventDefault();
            window.history.pushState({}, '', href);
            document.body.innerHTML = '<div style="text-align: center; padding: 3rem;"><h1>Loading...</h1><p>Please wait while we load the page.</p><p><a href="/">Return Home</a></p></div>';
            setTimeout(function() {
              window.location.href = href;
            }, 500);
          }
        });
      });
      
      window.addEventListener('popstate', function() {
        window.location.reload();
      });
    });
  </script>
</body>
</html>
EOL

# Create mixtures directory and index page
mkdir -p netlify-static-site/out/mixtures
cat > netlify-static-site/out/mixtures/index.html << EOL
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Mixtures | CryoProtect</title>
  <style>
    body {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans', 'Helvetica Neue', sans-serif;
      line-height: 1.6;
      color: #333;
      max-width: 800px;
      margin: 0 auto;
      padding: 2rem;
    }
    h1 {
      font-size: 2rem;
      margin-bottom: 1rem;
      color: #0070f3;
    }
    .card {
      border: 1px solid #eaeaea;
      border-radius: 0.5rem;
      padding: 1.5rem;
      margin: 1rem 0;
    }
    .button {
      display: inline-block;
      background-color: #0070f3;
      color: white;
      padding: 0.5rem 1rem;
      text-decoration: none;
      border-radius: 0.25rem;
      font-weight: 500;
      margin-right: 0.5rem;
    }
    .home-link {
      display: inline-flex;
      align-items: center;
      color: #0070f3;
      text-decoration: none;
      margin-bottom: 1rem;
    }
    .mixture-list {
      list-style: none;
      padding: 0;
    }
    .mixture-item {
      padding: 0.75rem;
      border-bottom: 1px solid #eaeaea;
    }
    .mixture-link {
      color: #0070f3;
      text-decoration: none;
      display: block;
    }
    .mixture-link:hover {
      text-decoration: underline;
    }
  </style>
</head>
<body>
  <a href="/" class="home-link">
    <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M19 12H5M12 19l-7-7 7-7"/></svg>
    <span style="margin-left: 0.5rem;">Home</span>
  </a>
  
  <h1>Mixtures</h1>
  <p>Browse and analyze cryoprotectant mixtures.</p>
  
  <div class="card">
    <h2>Featured Mixtures</h2>
    <ul class="mixture-list">
      <li class="mixture-item">
        <a href="/mixtures/1" class="mixture-link">DMSO/Glycerol (10%/10%)</a>
        <p>Standard cryopreservation mixture for cell storage</p>
      </li>
      <li class="mixture-item">
        <a href="/mixtures/2" class="mixture-link">DMSO/Trehalose (15%/5%)</a>
        <p>Enhanced cell viability for sensitive tissue samples</p>
      </li>
      <li class="mixture-item">
        <a href="/mixtures/3" class="mixture-link">Propylene Glycol/Ethylene Glycol (7.5%/7.5%)</a>
        <p>Optimized for embryo cryopreservation</p>
      </li>
    </ul>
    
    <div style="margin-top: 1.5rem;">
      <a href="/mixtures/new" class="button">Create New Mixture</a>
    </div>
  </div>
  
  <p style="margin-top: 2rem; color: #666; font-size: 0.875rem;">
    For the full mixture analysis experience, please use the main application with JavaScript enabled.
  </p>

  <script>
    // Client-side navigation handler
    document.addEventListener('DOMContentLoaded', function() {
      document.querySelectorAll('a').forEach(function(link) {
        link.addEventListener('click', function(e) {
          var href = this.getAttribute('href');
          
          if (href && href.indexOf('/') === 0 && href.indexOf('//') !== 0) {
            e.preventDefault();
            window.history.pushState({}, '', href);
            document.body.innerHTML = '<div style="text-align: center; padding: 3rem;"><h1>Loading...</h1><p>Please wait while we load the page.</p><p><a href="/">Return Home</a></p></div>';
            setTimeout(function() {
              window.location.href = href;
            }, 500);
          }
        });
      });
      
      window.addEventListener('popstate', function() {
        window.location.reload();
      });
    });
  </script>
</body>
</html>
EOL

# Create properties directory and index page
mkdir -p netlify-static-site/out/properties
cat > netlify-static-site/out/properties/index.html << EOL
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Properties | CryoProtect</title>
  <style>
    body {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans', 'Helvetica Neue', sans-serif;
      line-height: 1.6;
      color: #333;
      max-width: 800px;
      margin: 0 auto;
      padding: 2rem;
    }
    h1 {
      font-size: 2rem;
      margin-bottom: 1rem;
      color: #0070f3;
    }
    .card {
      border: 1px solid #eaeaea;
      border-radius: 0.5rem;
      padding: 1.5rem;
      margin: 1rem 0;
    }
    .button {
      display: inline-block;
      background-color: #0070f3;
      color: white;
      padding: 0.5rem 1rem;
      text-decoration: none;
      border-radius: 0.25rem;
      font-weight: 500;
      margin-right: 0.5rem;
    }
    .home-link {
      display: inline-flex;
      align-items: center;
      color: #0070f3;
      text-decoration: none;
      margin-bottom: 1rem;
    }
    .property-list {
      list-style: none;
      padding: 0;
    }
    .property-item {
      padding: 0.75rem;
      border-bottom: 1px solid #eaeaea;
    }
  </style>
</head>
<body>
  <a href="/" class="home-link">
    <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M19 12H5M12 19l-7-7 7-7"/></svg>
    <span style="margin-left: 0.5rem;">Home</span>
  </a>
  
  <h1>Properties</h1>
  <p>Explore molecular properties relevant to cryoprotection.</p>
  
  <div class="card">
    <h2>Featured Properties</h2>
    <ul class="property-list">
      <li class="property-item">
        <h3>Molecular Weight</h3>
        <p>The mass of a molecule, measured in atomic mass units (g/mol)</p>
      </li>
      <li class="property-item">
        <h3>LogP</h3>
        <p>The logarithm of the partition coefficient between octanol and water, indicating lipophilicity</p>
      </li>
      <li class="property-item">
        <h3>Glass Transition Temperature</h3>
        <p>Temperature at which a substance transitions from a crystalline to an amorphous state</p>
      </li>
      <li class="property-item">
        <h3>Hydrogen Bond Donors/Acceptors</h3>
        <p>Number of hydrogen bond donors and acceptors in a molecule</p>
      </li>
    </ul>
  </div>
  
  <div class="card">
    <h2>Property Search</h2>
    <p>Use our property explorer to find molecules with specific properties.</p>
    <p style="color: #666;">
      Property search requires JavaScript to be enabled. Please use the main application for the full experience.
    </p>
    <div style="margin-top: 1rem;">
      <a href="/molecules" class="button">Browse Molecules</a>
    </div>
  </div>
  
  <p style="margin-top: 2rem; color: #666; font-size: 0.875rem;">
    For the full property explorer experience, please use the main application with JavaScript enabled.
  </p>

  <script>
    // Client-side navigation handler
    document.addEventListener('DOMContentLoaded', function() {
      document.querySelectorAll('a').forEach(function(link) {
        link.addEventListener('click', function(e) {
          var href = this.getAttribute('href');
          
          if (href && href.indexOf('/') === 0 && href.indexOf('//') !== 0) {
            e.preventDefault();
            window.history.pushState({}, '', href);
            document.body.innerHTML = '<div style="text-align: center; padding: 3rem;"><h1>Loading...</h1><p>Please wait while we load the page.</p><p><a href="/">Return Home</a></p></div>';
            setTimeout(function() {
              window.location.href = href;
            }, 500);
          }
        });
      });
      
      window.addEventListener('popstate', function() {
        window.location.reload();
      });
    });
  </script>
</body>
</html>
EOL

# Create Netlify redirects file
echo "📝 Creating _redirects file for Netlify..."
cat > netlify-static-site/out/_redirects << EOL
# API routes redirect to Heroku backend
/api/*  https://cryoprotect-8030e4025428.herokuapp.com/api/:splat  200

# Handle dynamic routes for molecules
/molecules/[0-9]*  /molecules/detail.html  200
/molecules/*       /molecules/detail.html  200

# Handle dynamic routes for mixtures
/mixtures/[0-9]*  /mixtures/index.html  200
/mixtures/*       /mixtures/index.html  200

# SPA fallback for all other routes
/*  /index.html  200
EOL

# Copy the static files to the deployment directory
echo "📦 Copying static files to the deployment directory..."
rm -rf out
cp -r netlify-static-site/out ./out

echo "✅ Static site created successfully at ./out"
echo "🚀 Ready for deployment to Netlify!"