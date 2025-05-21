// Client-side verification for protection bypass
document.addEventListener('DOMContentLoaded', () => {
  const outputDiv = document.getElementById('output');
  const testBtn = document.getElementById('test-btn');
  
  function log(message, isError = false) {
    const entry = document.createElement('div');
    entry.textContent = message;
    if (isError) {
      entry.style.color = 'red';
    } else if (message.includes('SUCCESS')) {
      entry.style.color = 'green';
      entry.style.fontWeight = 'bold';
    }
    outputDiv.appendChild(entry);
  }

  function testEndpoint(url, withToken = false, useHeader = true) {
    const options = {
      method: 'GET',
      headers: {
        'Accept': 'application/json'
      }
    };
    
    let displayUrl = url;
    
    if (withToken) {
      const token = 'TAt23KbtFE8dkZobJU3hpgTP4L5ja07V';
      if (useHeader) {
        options.headers['x-protection-bypass'] = token;
        log(`Testing with token in header: ${url}`);
      } else {
        url = `${url}${url.includes('?') ? '&' : '?'}bypass=${token}`;
        displayUrl = url;
        log(`Testing with token in query param: ${url}`);
      }
    } else {
      log(`Testing without token: ${url}`);
    }
    
    return fetch(url, options)
      .then(response => {
        if (response.ok) {
          log(`✅ SUCCESS: ${response.status} - ${displayUrl}`);
          return response.text();
        } else {
          log(`❌ FAILURE: ${response.status} - ${displayUrl}`, true);
          return response.text().then(text => {
            log(`Response: ${text.substring(0, 100)}...`, true);
            throw new Error(`Request failed with status ${response.status}`);
          });
        }
      })
      .catch(error => {
        log(`ERROR: ${error.message}`, true);
      });
  }
  
  testBtn.addEventListener('click', () => {
    const baseUrl = document.getElementById('base-url').value;
    if (!baseUrl) {
      log('Please enter a base URL', true);
      return;
    }
    
    // Clear previous results
    outputDiv.innerHTML = '';
    
    const apiUrl = `${baseUrl}/api/v1/health`;
    
    // Test without token (should fail)
    testEndpoint(apiUrl)
      .finally(() => {
        // Test with token in header (should succeed)
        testEndpoint(apiUrl, true, true)
          .finally(() => {
            // Test with token in query param (should succeed)
            testEndpoint(apiUrl, true, false);
          });
      });
  });
});