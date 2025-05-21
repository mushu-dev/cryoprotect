#!/usr/bin/env python3
"""
Debug script to test all configured providers and API keys
"""

import os
from dotenv import load_dotenv
from os_computer_use import providers

# Load environment variables from .env file
load_dotenv()

def main():
    print("Testing all providers to identify issues:")
    
    # Test API keys
    api_keys = {
        "E2B_API_KEY": os.getenv("E2B_API_KEY"),
        "ANTHROPIC_API_KEY": os.getenv("ANTHROPIC_API_KEY"),
        "OPENAI_API_KEY": os.getenv("OPENAI_API_KEY"),
        "GROQ_API_KEY": os.getenv("GROQ_API_KEY"),
        "OPENROUTER_API_KEY": os.getenv("OPENROUTER_API_KEY"),
        "FIREWORKS_API_KEY": os.getenv("FIREWORKS_API_KEY"),
        "LLAMA_API_KEY": os.getenv("LLAMA_API_KEY"),
        "DEEPSEEK_API_KEY": os.getenv("DEEPSEEK_API_KEY"),
        "GEMINI_API_KEY": os.getenv("GEMINI_API_KEY"),
        "MISTRAL_API_KEY": os.getenv("MISTRAL_API_KEY"),
        "MOONSHOT_API_KEY": os.getenv("MOONSHOT_API_KEY"),
        "HF_TOKEN": os.getenv("HF_TOKEN")
    }
    
    print("\nAPI Keys:")
    for key, value in api_keys.items():
        status = "✅ Set" if value else "❌ Missing"
        print(f"  {key}: {status}")
    
    print("\nTesting OSAtlasProvider:")
    try:
        osatlas = providers.OSAtlasProvider()
        print("  OSAtlasProvider initialized without error")
    except Exception as e:
        print(f"  Error initializing OSAtlasProvider: {e}")
    
    print("\nTesting AnthropicProvider:")
    try:
        anthropic = providers.AnthropicProvider("claude-3.5-sonnet")
        print("  AnthropicProvider initialized without error")
    except Exception as e:
        print(f"  Error initializing AnthropicProvider: {e}")
    
    print("\nTesting OpenAIProvider:")
    try:
        openai = providers.OpenAIProvider("gpt-4o")
        print("  OpenAIProvider initialized without error")
    except Exception as e:
        print(f"  Error initializing OpenAIProvider: {e}")
    
    print("\nTesting GroqProvider:")
    try:
        groq = providers.GroqProvider("llama-3.3")
        print("  GroqProvider initialized without error")
    except Exception as e:
        print(f"  Error initializing GroqProvider: {e}")
    
    print("\nTesting OpenRouterProvider:")
    try:
        openrouter = providers.OpenRouterProvider("qwen-2.5-vl")
        print("  OpenRouterProvider initialized without error")
    except Exception as e:
        print(f"  Error initializing OpenRouterProvider: {e}")

if __name__ == "__main__":
    main()