# Using Anthropic for both vision and action models
from os_computer_use import providers

grounding_model = providers.OSAtlasProvider()
vision_model = providers.AnthropicProvider("claude-3.5-sonnet")
action_model = providers.AnthropicProvider("claude-3.5-sonnet")
