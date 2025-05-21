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
import type * as auth_clerk from "../auth/clerk.js";
import type * as auth_index from "../auth/index.js";
import type * as auth_users from "../auth/users.js";
import type * as dataSources_dataSources from "../dataSources/dataSources.js";
import type * as dataSources_helpers from "../dataSources/helpers.js";
import type * as dataSources_index from "../dataSources/index.js";
import type * as dataSources_types from "../dataSources/types.js";
import type * as dataSources_validation from "../dataSources/validation.js";
import type * as experiments_experimentResults from "../experiments/experimentResults.js";
import type * as experiments_experiments from "../experiments/experiments.js";
import type * as experiments_helpers from "../experiments/helpers.js";
import type * as experiments_index from "../experiments/index.js";
import type * as experiments_types from "../experiments/types.js";
import type * as experiments_validation from "../experiments/validation.js";
import type * as mixtures_helpers from "../mixtures/helpers.js";
import type * as mixtures_index from "../mixtures/index.js";
import type * as mixtures_mixtureComponents from "../mixtures/mixtureComponents.js";
import type * as mixtures_mixtures from "../mixtures/mixtures.js";
import type * as mixtures_types from "../mixtures/types.js";
import type * as mixtures_validation from "../mixtures/validation.js";
import type * as molecules_create from "../molecules/create.js";
import type * as molecules_delete from "../molecules/delete.js";
import type * as molecules_helpers from "../molecules/helpers.js";
import type * as molecules_index from "../molecules/index.js";
import type * as molecules_mutations from "../molecules/mutations.js";
import type * as molecules_queries from "../molecules/queries.js";
import type * as molecules_query from "../molecules/query.js";
import type * as molecules_types from "../molecules/types.js";
import type * as molecules_update from "../molecules/update.js";
import type * as molecules_validation from "../molecules/validation.js";
import type * as predictions_helpers from "../predictions/helpers.js";
import type * as predictions_index from "../predictions/index.js";
import type * as predictions_predictions from "../predictions/predictions.js";
import type * as predictions_scientificModels from "../predictions/scientificModels.js";
import type * as predictions_types from "../predictions/types.js";
import type * as predictions_validation from "../predictions/validation.js";
import type * as properties_helpers from "../properties/helpers.js";
import type * as properties_molecularProperties from "../properties/molecularProperties.js";
import type * as properties_propertyTypes from "../properties/propertyTypes.js";
import type * as properties_types from "../properties/types.js";
import type * as properties_validation from "../properties/validation.js";
import type * as references_crossReferences from "../references/crossReferences.js";
import type * as references_helpers from "../references/helpers.js";
import type * as references_index from "../references/index.js";
import type * as references_synonyms from "../references/synonyms.js";
import type * as references_types from "../references/types.js";
import type * as references_validation from "../references/validation.js";
import type * as schema_convex_schema from "../schema/convex_schema.js";
import type * as users_helpers from "../users/helpers.js";
import type * as users_index from "../users/index.js";
import type * as users_projects from "../users/projects.js";
import type * as users_teamMembers from "../users/teamMembers.js";
import type * as users_types from "../users/types.js";
import type * as users_users from "../users/users.js";
import type * as users_validation from "../users/validation.js";
import type * as utils_common from "../utils/common.js";

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
  "auth/clerk": typeof auth_clerk;
  "auth/index": typeof auth_index;
  "auth/users": typeof auth_users;
  "dataSources/dataSources": typeof dataSources_dataSources;
  "dataSources/helpers": typeof dataSources_helpers;
  "dataSources/index": typeof dataSources_index;
  "dataSources/types": typeof dataSources_types;
  "dataSources/validation": typeof dataSources_validation;
  "experiments/experimentResults": typeof experiments_experimentResults;
  "experiments/experiments": typeof experiments_experiments;
  "experiments/helpers": typeof experiments_helpers;
  "experiments/index": typeof experiments_index;
  "experiments/types": typeof experiments_types;
  "experiments/validation": typeof experiments_validation;
  "mixtures/helpers": typeof mixtures_helpers;
  "mixtures/index": typeof mixtures_index;
  "mixtures/mixtureComponents": typeof mixtures_mixtureComponents;
  "mixtures/mixtures": typeof mixtures_mixtures;
  "mixtures/types": typeof mixtures_types;
  "mixtures/validation": typeof mixtures_validation;
  "molecules/create": typeof molecules_create;
  "molecules/delete": typeof molecules_delete;
  "molecules/helpers": typeof molecules_helpers;
  "molecules/index": typeof molecules_index;
  "molecules/mutations": typeof molecules_mutations;
  "molecules/queries": typeof molecules_queries;
  "molecules/query": typeof molecules_query;
  "molecules/types": typeof molecules_types;
  "molecules/update": typeof molecules_update;
  "molecules/validation": typeof molecules_validation;
  "predictions/helpers": typeof predictions_helpers;
  "predictions/index": typeof predictions_index;
  "predictions/predictions": typeof predictions_predictions;
  "predictions/scientificModels": typeof predictions_scientificModels;
  "predictions/types": typeof predictions_types;
  "predictions/validation": typeof predictions_validation;
  "properties/helpers": typeof properties_helpers;
  "properties/molecularProperties": typeof properties_molecularProperties;
  "properties/propertyTypes": typeof properties_propertyTypes;
  "properties/types": typeof properties_types;
  "properties/validation": typeof properties_validation;
  "references/crossReferences": typeof references_crossReferences;
  "references/helpers": typeof references_helpers;
  "references/index": typeof references_index;
  "references/synonyms": typeof references_synonyms;
  "references/types": typeof references_types;
  "references/validation": typeof references_validation;
  "schema/convex_schema": typeof schema_convex_schema;
  "users/helpers": typeof users_helpers;
  "users/index": typeof users_index;
  "users/projects": typeof users_projects;
  "users/teamMembers": typeof users_teamMembers;
  "users/types": typeof users_types;
  "users/users": typeof users_users;
  "users/validation": typeof users_validation;
  "utils/common": typeof utils_common;
}>;
export declare const api: FilterApi<
  typeof fullApi,
  FunctionReference<any, "public">
>;
export declare const internal: FilterApi<
  typeof fullApi,
  FunctionReference<any, "internal">
>;
