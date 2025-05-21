/* eslint-disable */
/**
 * Generated `api` utility.
 *
 * THIS CODE IS AUTOMATICALLY GENERATED.
 *
 * To regenerate, run `npx convex dev`.
 * @module
 */

import type {
  ApiFromModules,
  FilterApi,
  FunctionReference,
} from "convex/server";
import type * as api_auth_signin from "../api/auth/signin.js";
import type * as api_auth_signout from "../api/auth/signout.js";
import type * as api_auth_signup from "../api/auth/signup.js";
import type * as api_delete from "../api/delete.js";
import type * as api_insert from "../api/insert.js";
import type * as api_query from "../api/query.js";
import type * as api_update from "../api/update.js";
import type * as http_api from "../http/api.js";

/**
 * A utility for referencing Convex functions in your app's API.
 *
 * Usage:
 * ```js
 * const myFunctionReference = api.myModule.myFunction;
 * ```
 */
declare const fullApi: ApiFromModules<{
  "api/auth/signin": typeof api_auth_signin;
  "api/auth/signout": typeof api_auth_signout;
  "api/auth/signup": typeof api_auth_signup;
  "api/delete": typeof api_delete;
  "api/insert": typeof api_insert;
  "api/query": typeof api_query;
  "api/update": typeof api_update;
  "http/api": typeof http_api;
}>;
export declare const api: FilterApi<
  typeof fullApi,
  FunctionReference<any, "public">
>;
export declare const internal: FilterApi<
  typeof fullApi,
  FunctionReference<any, "internal">
>;
