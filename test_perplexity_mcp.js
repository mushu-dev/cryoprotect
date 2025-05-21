// Simple script to test Perplexity MCP server
const { spawnSync } = require('child_process');
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

// Try to run the server process
console.log(`Starting server with command: ${serverConfig.command} ${serverConfig.args.join(' ')}`);

// Prepare environment variables
const env = { ...process.env };
if (serverConfig.env) {
  Object.assign(env, serverConfig.env);
}

// Run the server process
try {
  const result = spawnSync(serverConfig.command, serverConfig.args, {
    env,
    stdio: 'pipe',
    encoding: 'utf8',
    timeout: 5000 // 5 seconds timeout
  });

  if (result.error) {
    console.error('Error starting server:', result.error);
    process.exit(1);
  }

  console.log('Server process started successfully');
  console.log('Stdout:', result.stdout);
  console.log('Stderr:', result.stderr);
  
  if (result.status !== 0) {
    console.error(`Server process exited with code ${result.status}`);
    process.exit(1);
  }
  
  console.log('Server test completed successfully');
} catch (error) {
  console.error('Failed to start server process:', error);
  process.exit(1);
}