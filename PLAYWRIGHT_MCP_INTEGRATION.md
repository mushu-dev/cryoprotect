# Playwright MCP Integration

This directory contains the integration between the CryoProtect project and the Playwright MCP (Model Control Protocol) service. This integration allows us to use Playwright for browser automation through Claude Code's MCP interface.

## Components

1. **mcp-playwright-final.sh**: The main script that handles commands from Claude Code MCP and executes them in a containerized Playwright environment.

2. **test-cryoprotect.js**: A sample Playwright test script that demonstrates how to test the CryoProtect application.

3. **playwright.config.js**: Configuration file for Playwright tests.

4. **run-playwright-test.sh**: A script to run Playwright tests in a container environment without requiring MCP.

## Usage

### Through Claude Code MCP

When using Claude Code, the MCP Playwright integration is used automatically when you call the appropriate MCP commands:

```
mcp__playwright-docker__browser_navigate https://cryoprotect.app
mcp__playwright-docker__browser_take_screenshot https://cryoprotect.app screenshot.png
mcp__playwright-docker__browser_snapshot https://cryoprotect.app
```

### Manually Testing with the Script

You can also use the scripts directly:

```bash
# Navigate to a URL and print information
./mcp-playwright-final.sh browser_navigate https://cryoprotect.app

# Take a screenshot of a website
./mcp-playwright-final.sh browser_take_screenshot https://cryoprotect.app screenshot.png

# Get accessibility snapshot of a website
./mcp-playwright-final.sh browser_snapshot https://cryoprotect.app

# Check if the integration is working
./mcp-playwright-final.sh status
```

### Running the Test Script

To run the sample test script:

```bash
./run-playwright-test.sh
```

This will execute the test in a container environment and generate a report.

## Benefits

- **Cross-platform compatibility**: Works on any system that supports Docker or Podman
- **No local dependencies**: All required libraries are contained within the container
- **Isolated environment**: Tests run in a clean, controlled environment
- **Seamless MCP integration**: Works transparently with Claude Code's MCP interface
- **Browser automation**: Supports navigation, screenshots, accessibility snapshots, and more

## Troubleshooting

If you encounter any issues:

1. Make sure Docker or Podman is installed and running on your system
2. Check that the container has internet access
3. Verify that the scripts have executable permissions (`chmod +x *.sh`)
4. Use the `status` command to verify the integration is working correctly

For more detailed information, see the `CLAUDE.md` file in the root directory.