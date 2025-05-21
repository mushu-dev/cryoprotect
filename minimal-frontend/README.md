# CryoProtect Minimal Frontend

This is a minimal frontend implementation for the CryoProtect application. It provides a basic structure and placeholder pages for the main sections of the application.

## Features

- Responsive design for all device sizes
- Basic navigation structure
- Placeholder pages for main sections:
  - Molecules
  - Mixtures
  - Protocols
  - Experiments
- About page with project information
- Custom 404 error page
- Consistent layout and styling across all pages

## Getting Started

### Prerequisites

- Node.js 14.x or later
- npm 6.x or later

### Installation

1. Clone the repository:
   ```
   git clone https://github.com/your-org/cryoprotect.git
   cd cryoprotect/minimal-frontend
   ```

2. Install dependencies:
   ```
   npm install
   ```

3. Run the development server:
   ```
   npm run dev
   ```

4. Open [http://localhost:3000](http://localhost:3000) with your browser to see the result.

## Deployment

This project is configured for deployment to Netlify. To deploy:

1. Build the project:
   ```
   npm run build
   npm run export
   ```

2. Deploy to Netlify:
   ```
   ./deploy-to-netlify.sh
   ```

Alternatively, you can set up automatic deployments by connecting your GitHub repository to Netlify.

## Project Structure

- `pages/`: Contains all page components
- `components/`: Reusable UI components
- `styles/`: CSS styles
- `public/`: Static assets

## Next Steps

- Connect to the backend API
- Implement actual data fetching and display
- Add user authentication
- Implement full protocol and experiment management features
- Add molecular visualization tools

## License

This project is licensed under the MIT License - see the LICENSE file for details.