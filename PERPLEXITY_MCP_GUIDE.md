# Perplexity MCP Server Guide

This guide explains how to use the Perplexity MCP server with Claude Code. The server provides powerful research and reasoning capabilities powered by Perplexity's specialized AI models.

## Setup and Configuration

### 1. Prerequisites
- Node.js installed (v18+)
- Perplexity API key from [perplexity.ai/settings/api](https://www.perplexity.ai/settings/api)
- The perplexity-mcp repository cloned locally

### 2. Building the Server
```bash
cd /home/mushu/Projects/cryoprotect/perplexity-mcp
npm run build
```

### 3. Claude Code MCP Configuration
The MCP configuration file is located at `~/.config/mcp/claude_desktop_config.json`. Ensure it contains the following:

```json
{
  "mcpServers": {
    "perplexity": {
      "command": "node",
      "args": ["/home/mushu/Projects/cryoprotect/perplexity-mcp/build/index.js"],
      "env": {
        "PERPLEXITY_API_KEY": "your-api-key-here"
      },
      "disabled": false
    }
  }
}
```

Replace `your-api-key-here` with your actual Perplexity API key.

## Available Tools

The Perplexity MCP server provides three powerful tools with automatic query classification:

### 1. Search (Sonar Pro)
Best for simple queries and basic information lookup. Use this for straightforward questions that need concise answers.

Example:
```
What is the capital of France?
```

### 2. Reason (Sonar Reasoning Pro)
Handles complex tasks requiring deeper analysis. Perfect for explanations, comparisons, and problem-solving.

Example:
```
Compare and contrast REST and GraphQL APIs, explaining their pros and cons.
```

### 3. Deep Research (Sonar Deep Research)
Conducts comprehensive research on complex topics. Ideal for in-depth analysis requiring detailed investigation.

Example:
```
What are the latest advancements in quantum computing for cryptography, focusing on post-quantum algorithms?
```

## Intelligent Model Selection

The server automatically analyzes query complexity to route requests to the most appropriate model:

1. Simple queries → Sonar Pro (basic information lookup)
2. Complex queries → Sonar Reasoning Pro (how/why questions, comparisons)
3. Research queries → Sonar Deep Research (comprehensive analysis, detailed investigations)

You can override the automatic selection by using the `force_model` parameter.

## Testing the Server

Use the provided test scripts to verify the server is working correctly:

```bash
# Test server startup
node test_perplexity_mcp.js

# Test tool listing
node test_claude_perplexity.js

# Test calling a specific tool
node test_perplexity_tool.js
```

## Usage in Claude Code

Once configured, the Perplexity MCP server will be available in Claude Code through the MCP tools interface. You can use it to:

1. Research complex topics
2. Get detailed explanations of concepts
3. Compare and analyze different approaches
4. Answer factual questions with the latest information

## Troubleshooting

If you encounter issues with the Perplexity MCP server:

1. Check that the API key is valid and correctly configured
2. Ensure the server is built successfully (`npm run build`)
3. Verify the configuration file path is correct
4. Check for any error messages in the server output
5. Try restarting Claude Code after making configuration changes

For more detailed information, refer to the [official repository](https://github.com/DaInfernalCoder/perplexity-mcp).