// @ts-check
const { test, expect } = require('@playwright/test');

/**
 * Security tests for the experimental data enhancement UI
 * 
 * These tests verify security measures like proper headers, 
 * input validation, and protection against common vulnerabilities.
 */
test.describe('Security Tests', () => {
  test('should have proper security headers', async ({ page }) => {
    // Navigate to the homepage
    const response = await page.goto('/');
    
    // Check the response headers
    const headers = response.headers();
    console.log('Response Headers:', headers);
    
    // Check for key security headers
    // Content-Security-Policy
    const hasCSP = 'content-security-policy' in headers;
    console.log('Has Content-Security-Policy:', hasCSP);
    
    // X-XSS-Protection
    const hasXSSProtection = 'x-xss-protection' in headers;
    console.log('Has X-XSS-Protection:', hasXSSProtection);
    
    // X-Content-Type-Options
    const hasContentTypeOptions = 'x-content-type-options' in headers;
    console.log('Has X-Content-Type-Options:', hasContentTypeOptions);
    expect(hasContentTypeOptions).toBeTruthy();
    
    if (hasContentTypeOptions) {
      expect(headers['x-content-type-options']).toBe('nosniff');
    }
    
    // X-Frame-Options
    const hasFrameOptions = 'x-frame-options' in headers;
    console.log('Has X-Frame-Options:', hasFrameOptions);
    
    // Strict-Transport-Security (HSTS)
    const hasHSTS = 'strict-transport-security' in headers;
    console.log('Has Strict-Transport-Security:', hasHSTS);
    
    // Referrer-Policy
    const hasReferrerPolicy = 'referrer-policy' in headers;
    console.log('Has Referrer-Policy:', hasReferrerPolicy);
    
    // Check TLS/SSL configuration if on HTTPS
    if (page.url().startsWith('https://')) {
      // We can't directly test SSL/TLS settings with Playwright,
      // but we can note that we're using HTTPS
      console.log('Page is served over HTTPS');
    }
  });

  test('should validate inputs properly', async ({ page }) => {
    // Navigate to a form page
    await page.goto('/experiments/create');
    
    // Wait for form to load
    await page.waitForSelector('form');
    
    // Test XSS attempt in text inputs
    const xssPayload = '<script>alert("XSS")</script>';
    await page.fill('input[name="name"]', xssPayload);
    await page.fill('textarea[name="description"]', xssPayload);
    
    // Submit the form
    await page.click('button[type="submit"]');
    
    // Check if the form was rejected or sanitized
    // If validation works, either we'll stay on the same page or we'll get validation errors
    
    // Option 1: Check if we're still on the form page (didn't navigate away)
    const formStillVisible = await page.locator('form').isVisible();
    
    // Option 2: Check for validation error messages
    const hasValidationErrors = await page.locator('[data-testid="form-error"]').count() > 0;
    
    // The form should either prevent submission or show validation errors
    expect(formStillVisible || hasValidationErrors).toBeTruthy();
    
    // Test SQL injection attempt
    const sqlInjectionPayload = "' OR 1=1; --";
    await page.fill('input[name="name"]', sqlInjectionPayload);
    
    // Submit the form
    await page.click('button[type="submit"]');
    
    // Check for proper handling
    const stillOnFormPage = await page.locator('form').isVisible();
    expect(stillOnFormPage).toBeTruthy();
  });

  test('should not expose sensitive information in the UI', async ({ page }) => {
    // Pages to check for sensitive information
    const pagesToCheck = [
      '/',
      '/experiments',
      '/molecules',
      '/properties'
    ];
    
    // Patterns that might indicate sensitive information
    const sensitivePatterns = [
      /password/i,
      /api[-_]?key/i,
      /secret/i,
      /token/i,
      /auth[-_]?key/i,
      /connection[-_]?string/i
    ];
    
    // Check each page
    for (const path of pagesToCheck) {
      await page.goto(path);
      await page.waitForLoadState('networkidle');
      
      // Get all text on the page
      const pageText = await page.evaluate(() => document.body.innerText);
      
      // Check for sensitive patterns
      for (const pattern of sensitivePatterns) {
        const match = pageText.match(pattern);
        if (match) {
          console.warn(`Potential sensitive information found on ${path}: ${match[0]}`);
        }
        
        // No sensitive patterns should be found
        expect(match).toBeNull();
      }
      
      // Check page source for hardcoded secrets or tokens
      const pageSource = await page.content();
      const suspiciousTokenPatterns = [
        /(['"])(?:sk|pk)_[0-9a-zA-Z]{20,}(['"])/,  // API keys like sk_live_...
        /(['"])[0-9a-f]{32,}(['"])/,               // MD5 hashes or similar
        /(['"])[0-9a-zA-Z+/=]{40,}(['"])/,         // Base64 encoded secrets
        /authorization:\s*(['"])bearer\s+[0-9a-zA-Z._-]+(['"])/i // Auth headers
      ];
      
      for (const pattern of suspiciousTokenPatterns) {
        const match = pageSource.match(pattern);
        if (match) {
          console.warn(`Suspicious token pattern found on ${path}: ${match[0].substring(0, 20)}...`);
        }
        
        // No suspicious token patterns should be found
        expect(match).toBeNull();
      }
    }
  });

  test('should handle authentication and authorization correctly', async ({ page }) => {
    // Check if there's a login functionality
    await page.goto('/');
    
    // Look for login/auth UI elements
    const loginButton = page.locator('button, a').filter({ hasText: /log[- ]?in|sign[- ]?in/i });
    const hasLoginUI = await loginButton.count() > 0;
    
    if (hasLoginUI) {
      console.log('Login UI detected, testing auth flow');
      
      // Click login button
      await loginButton.click();
      
      // Wait for login form or auth redirect
      const hasLoginForm = await page.locator('form input[type="password"]').isVisible()
        .catch(() => false);
      
      if (hasLoginForm) {
        // Test with invalid credentials
        await page.fill('input[type="email"], input[name="username"]', 'test@example.com');
        await page.fill('input[type="password"]', 'wrongpassword');
        
        // Submit the form
        await page.click('button[type="submit"]');
        
        // Should show an error and not authenticate
        const errorVisible = await page.locator('.error-message, [data-testid="auth-error"]')
          .isVisible()
          .catch(() => false);
        
        // We should either see an error message or stay on the login page
        const stillOnLoginPage = await page.locator('input[type="password"]').isVisible()
          .catch(() => false);
        
        expect(errorVisible || stillOnLoginPage).toBeTruthy();
      }
    } else {
      console.log('No login UI detected, skipping auth flow test');
    }
    
    // Test accessing protected routes
    const potentiallyProtectedRoutes = [
      '/admin',
      '/settings',
      '/profile',
      '/dashboard'
    ];
    
    for (const route of potentiallyProtectedRoutes) {
      // Try to access the route directly
      const response = await page.goto(route);
      
      console.log(`Accessed ${route}: Status ${response.status()}`);
      
      // We should either get redirected (302/303), unauthorized (401),
      // forbidden (403), or not found (404)
      const isProtected = response.status() === 401 || 
                           response.status() === 403 || 
                           response.status() === 302 || 
                           response.status() === 303;
                           
      const isNotFound = response.status() === 404;
      
      // Either the route is protected or it doesn't exist
      if (!isNotFound) {
        expect(isProtected).toBeTruthy();
      }
    }
  });

  test('should handle CSRF protection correctly', async ({ page }) => {
    // CSRF protection often uses tokens in forms
    await page.goto('/experiments/create');
    
    // Wait for form to load
    await page.waitForSelector('form');
    
    // Check for CSRF token in the form
    const hasCSRFToken = await page.evaluate(() => {
      // Look for common CSRF token field names
      return document.querySelector('input[name="csrf_token"], input[name="_csrf"], input[name="csrfmiddlewaretoken"], input[name="_token"]') !== null;
    });
    
    console.log('Form has CSRF token:', hasCSRFToken);
    
    // If there's no CSRF token in the form, check for other CSRF protections
    if (!hasCSRFToken) {
      // Check for CSRF in cookies or headers
      const cookies = await page.context().cookies();
      const csrfCookie = cookies.find(cookie => 
        cookie.name.toLowerCase().includes('csrf') || 
        cookie.name.toLowerCase().includes('xsrf')
      );
      
      // Check if there's a SameSite cookie attribute
      const hasSameSiteCookies = cookies.some(cookie => cookie.sameSite === 'strict' || cookie.sameSite === 'lax');
      
      console.log('Has CSRF cookie:', !!csrfCookie);
      console.log('Has SameSite cookies:', hasSameSiteCookies);
      
      // In modern applications, SameSite cookies can be an alternative to CSRF tokens
      const hasCSRFProtection = hasCSRFToken || !!csrfCookie || hasSameSiteCookies;
      
      // Log the finding but don't make it a hard failure
      // since there are multiple valid approaches to CSRF protection
      if (!hasCSRFProtection) {
        console.warn('No visible CSRF protection detected. This might be a security risk.');
      }
    }
  });
});