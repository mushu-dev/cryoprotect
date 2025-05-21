import * as THREE from 'three';

declare module '@react-three/fiber' {
  interface ThreeElements {
    mesh: JSX.IntrinsicElements['mesh'] & { geometry?: THREE.BufferGeometry };
    bufferGeometry: JSX.IntrinsicElements['bufferGeometry'];
    lineSegments: JSX.IntrinsicElements['lineSegments'];
  }
}