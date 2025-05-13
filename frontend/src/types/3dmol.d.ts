declare module '3dmol' {
  export function createViewer(element: HTMLElement, config: {
    backgroundColor?: string;
    antialias?: boolean;
    width?: string | number;
    height?: string | number;
    [key: string]: any;
  }): any;
  
  export function generate3DStructure(model: any): void;
  
  export enum SpriteAlignment {
    topLeft,
    topCenter,
    topRight,
    centerLeft,
    center,
    centerRight,
    bottomLeft,
    bottomCenter,
    bottomRight
  }
}