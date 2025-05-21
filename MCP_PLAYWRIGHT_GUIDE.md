# MCP Playwright Guide

This guide explains how to use MCP Playwright on Linux systems, particularly Fedora, where there are library compatibility challenges.

## Solution Overview

We've implemented a robust containerized solution for MCP Playwright that works on any Linux system:

1. **Container-based approach**: Uses the official Microsoft Playwright Docker container
2. **Zero system modifications**: No need to install specific library versions or modify system libraries
3. **Complete compatibility**: Works with all MCP Playwright commands
4. **Easy to use**: Simple shell scripts to interact with the container
5. **JSON output**: Standardized JSON responses for easy parsing and integration
6. **Automatic setup**: Container is pulled and configured automatically as needed
7. **Multiple implementation options**: Choose the approach that best fits your needs

## How to Use

### Primary Solution: Containerized Playwright

The `mcp-playwright-final.sh` script provides a complete solution with these commands:

```bash
# Check if Playwright container is working properly
./mcp-playwright-final.sh status

# Navigate to a URL and get information
./mcp-playwright-final.sh browser_navigate https://example.com

# Take a screenshot of a website
./mcp-playwright-final.sh browser_take_screenshot https://example.com output.png

# Get accessibility snapshot of a website
./mcp-playwright-final.sh browser_snapshot https://example.com
```

### Testing the Solution

We provide test scripts to verify the solution works correctly:

```bash
# Test the containerized solution (recommended)
node test-final-playwright.js

# Test the direct container approach
node test-direct-playwright.js

# Test the original symbolic link approach
./mcp-cli.sh node test_mcp_services.js
```

## How It Works

Our solution addresses the Playwright compatibility issues that arise on newer Linux distributions:

1. **Problem**: Playwright expects specific library versions (libicu.so.66, libjpeg.so.8, etc.) but Fedora 42+ has newer versions
2. **Container Solution**: Use the Microsoft Playwright container which has all required libraries pre-installed
3. **Integration**: Shell scripts create temporary Node.js scripts that run inside the container
4. **File Sharing**: Directories are mounted to share files between host and container
5. **JSON Output**: All commands return standardized JSON responses for consistent parsing

## Alternative Approaches

We provide multiple solutions to ensure compatibility across different environments:

1. **mcp-playwright-final.sh**: Production-ready solution using containers (recommended)
2. **mcp-playwright-direct.sh**: Alternative implementation with direct container access
3. **mcp-cli.sh + library symlinks**: For systems without container support

## Requirements

- Podman or Docker installed on your system
- Internet access (first run will download the container image)
- Basic shell scripting knowledge for customization

## Troubleshooting

### Container Issues

- **Container won't start**: Check if podman/docker is running with `systemctl status podman`
- **Permission errors**: Run with sudo or add your user to the appropriate group
- **Network issues**: Check if container can access the internet

### MCP Integration Issues

- **MCP commands not working**: Verify MCP status with `mcp` in Claude Code CLI
- **Invalid JSON**: Check that commands are returning valid JSON output
- **File path issues**: Use absolute paths when specifying file locations

## Advanced Usage

### Adding Custom Playwright Commands

To add a new command to the container script:

1. Create a new function in `mcp-playwright-final.sh`
2. Create a JavaScript template that uses Playwright
3. Add the command to the main() function's case statement

Example for a new command:

```bash
browser_custom_command() {
  local url="$1"
  log "INFO" "Running custom command on URL: $url"
  
  # Create script
  cat > "$TEMP_DIR/custom.js" << EOF
const { chromium } = require('playwright');
(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  try {
    await page.goto('$url');
    // Custom logic here
    console.log(JSON.stringify({ success: true, data: { result: "custom data" } }));
  } catch (error) {
    console.log(JSON.stringify({ success: false, error: error.message }));
  } finally {
    await browser.close();
  }
})();
EOF

  # Run in container
  ensure_container
  run_in_container "cd /input && npm init -y && npm install playwright && node custom.js"
}
```

## Further Resources

- [Microsoft Playwright Documentation](https://playwright.dev/)
- [Playwright Docker Images](https://hub.docker.com/_/microsoft-playwright)
- [Podman Documentation](https://podman.io/docs)

## Support

For issues or questions about this solution, please create a GitHub issue or contact the development team.

---

This solution was developed to ensure full MCP Playwright functionality on Fedora and other modern Linux distributions.

✅ Tested on Fedora 42 with Claude Code