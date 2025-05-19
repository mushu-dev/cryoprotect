'use client'

import { SVGProps } from 'react'

// Flask icon
export function Flask(props: SVGProps<SVGSVGElement>) {
  return (
    <svg
      xmlns="http://www.w3.org/2000/svg"
      width="24"
      height="24"
      viewBox="0 0 24 24"
      fill="none"
      stroke="currentColor"
      strokeWidth="2"
      strokeLinecap="round"
      strokeLinejoin="round"
      {...props}
    >
      <path d="M9 3h6v3H9z"/>
      <path d="M10 9v14a1 1 0 0 0 1 1h2a1 1 0 0 0 1-1V9"/>
      <path d="M10 3v2"/>
      <path d="M14 3v2"/>
      <path d="M8.5 7h7"/>
    </svg>
  )
}

// FlaskConical icon
export function FlaskConical(props: SVGProps<SVGSVGElement>) {
  return (
    <svg
      xmlns="http://www.w3.org/2000/svg"
      width="24"
      height="24"
      viewBox="0 0 24 24"
      fill="none"
      stroke="currentColor"
      strokeWidth="2"
      strokeLinecap="round"
      strokeLinejoin="round"
      {...props}
    >
      <path d="M10 2v6l-2 3v9a2 2 0 0 0 2 2h4a2 2 0 0 0 2-2v-9l-2-3V2"/>
      <path d="M7 16h10"/>
    </svg>
  )
}

// ChartBar icon
export function ChartBar(props: SVGProps<SVGSVGElement>) {
  return (
    <svg
      xmlns="http://www.w3.org/2000/svg"
      width="24"
      height="24"
      viewBox="0 0 24 24"
      fill="none"
      stroke="currentColor"
      strokeWidth="2"
      strokeLinecap="round"
      strokeLinejoin="round"
      {...props}
    >
      <line x1="12" y1="20" x2="12" y2="10"/>
      <line x1="18" y1="20" x2="18" y2="4"/>
      <line x1="6" y1="20" x2="6" y2="16"/>
      <line x1="2" y1="20" x2="22" y2="20"/>
    </svg>
  )
}