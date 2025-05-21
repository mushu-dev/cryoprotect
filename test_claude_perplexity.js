// This script simulates what Claude Code would do to test an MCP server
const { spawn } = require('child_process');
const fs = require('fs');
const path = require('path');

// Get the configuration file path
const configPath = path.join(process.env.HOME, '.config', 'mcp', 'claude_desktop_config.json');
console.log(`Looking for config at: ${configPath}`);

// Check if the config file exists
if (!fs.existsSync(configPath)) {
  console.error('Config file not found:', configPath);
  process.exit(1);
}

// Read the config file
const config = JSON.parse(fs.readFileSync(configPath, 'utf8'));
console.log('Config loaded successfully');

// Check if perplexity server is configured
if (!config.mcpServers || !config.mcpServers.perplexity) {
  console.error('Perplexity server not configured in MCP config');
  process.exit(1);
}

const serverConfig = config.mcpServers.perplexity;
console.log('Perplexity server config found:', JSON.stringify(serverConfig, null, 2));

// Prepare environment variables
const env = { ...process.env };
if (serverConfig.env) {
  Object.assign(env, serverConfig.env);
}

// Start the server process
const serverProcess = spawn(serverConfig.command, serverConfig.args, {
  env,
  stdio: ['pipe', 'pipe', process.stderr]
});

// Prepare sample request for listing tools
const listToolsRequest = {
  jsonrpc: '2.0',
  id: 1,
  method: 'tools/list',
  params: {}
};

// Write the request to the server's stdin
serverProcess.stdin.write(JSON.stringify(listToolsRequest) + '\n');

// Set up timeout to kill the process after a while
const timeout = setTimeout(() => {
  console.error('Timeout waiting for response from server');
  serverProcess.kill();
  process.exit(1);
}, 10000);

// Collect response data
let responseData = '';
serverProcess.stdout.on('data', (data) => {
  responseData += data.toString();
  
  // Check if we have a complete response
  try {
    // Some MCP servers might output multiple lines, we need to find a valid JSON response
    const lines = responseData.split('\n').filter(line => line.trim());
    
    for (const line of lines) {
      try {
        const response = JSON.parse(line);
        
        if (response.id === 1) {
          // We got a valid response!
          clearTimeout(timeout);
          console.log('Received response from MCP server:');
          console.log(JSON.stringify(response, null, 2));
          
          // Check if the server returned tools
          if (response.result && response.result.tools) {
            console.log('Available tools:');
            response.result.tools.forEach((tool, index) => {
              console.log(`${index + 1}. ${tool.name}: ${tool.description}`);
            });
          }
          
          // Clean up and exit
          serverProcess.stdin.end();
          setTimeout(() => {
            serverProcess.kill();
            process.exit(0);
          }, 500);
          
          return;
        }
      } catch (e) {
        // Not a valid JSON line, continue
      }
    }
  } catch (e) {
    // Not a complete response yet
  }
});

serverProcess.on('error', (error) => {
  console.error('Error starting server:', error);
  clearTimeout(timeout);
  process.exit(1);
});

serverProcess.on('exit', (code, signal) => {
  if (code !== 0 && !signal) {
    console.error(`Server process exited with code ${code}`);
    clearTimeout(timeout);
    process.exit(1);
  }
});