# CryoProtect v2 Architecture Diagrams

## System Architecture Overview

```mermaid
graph TB
    subgraph "User Interface"
        UI[Web Browser]
    end

    subgraph "Application Layer"
        API[CryoProtect API]
        Auth[Authentication]
        RBAC[Role-Based Access Control]
    end

    subgraph "Data Layer"
        DB[(Supabase Database)]
        Cache[(Redis Cache)]
    end

    subgraph "Monitoring & Observability"
        Prom[Prometheus]
        Graf[Grafana]
        Alert[Alertmanager]
        NodeExp[Node Exporter]
    end

    subgraph "Logging"
        ES[Elasticsearch]
        LS[Logstash]
        KB[Kibana]
        FB[Filebeat]
    end

    subgraph "Backup"
        BS[Backup Scheduler]
        BM[Backup Manager]
        BR[Backup Repository]
    end

    UI --> API
    API --> Auth
    API --> RBAC
    API --> DB
    API --> Cache
    
    API --> Prom
    Prom --> Graf
    Prom --> Alert
    NodeExp --> Prom
    
    API --> LS
    LS --> ES
    ES --> KB
    FB --> LS
    
    BM --> DB
    BS --> BM
    BM --> BR
```

## Monitoring & Observability Architecture

```mermaid
graph TB
    subgraph "Application"
        API[CryoProtect API]
        Metrics[Prometheus Metrics]
        Middleware[Observability Middleware]
    end

    subgraph "Monitoring Stack"
        Prom[Prometheus]
        Graf[Grafana]
        Alert[Alertmanager]
        NodeExp[Node Exporter]
    end

    subgraph "Alerting"
        Email[Email Notifications]
        Slack[Slack Notifications]
    end

    subgraph "Dashboards"
        APIDash[API Performance]
        DBDash[Database Performance]
        SysDash[System Resources]
    end

    API --> Metrics
    API --> Middleware
    Metrics --> Prom
    NodeExp --> Prom
    Prom --> Graf
    Prom --> Alert
    Alert --> Email
    Alert --> Slack
    Graf --> APIDash
    Graf --> DBDash
    Graf --> SysDash
```

## Logging Architecture

```mermaid
graph TB
    subgraph "Application"
        API[CryoProtect API]
        Logger[Enhanced Logger]
        LogFiles[Log Files]
    end

    subgraph "ELK Stack"
        ES[Elasticsearch]
        LS[Logstash]
        KB[Kibana]
        FB[Filebeat]
    end

    subgraph "Log Analysis"
        Search[Search & Query]
        Viz[Visualizations]
        Dash[Dashboards]
    end

    API --> Logger
    Logger --> LogFiles
    LogFiles --> FB
    FB --> LS
    LS --> ES
    ES --> KB
    KB --> Search
    KB --> Viz
    KB --> Dash
```

## Backup Architecture

```mermaid
graph TB
    subgraph "Application"
        API[CryoProtect API]
        DB[(Supabase Database)]
    end

    subgraph "Backup System"
        BS[Backup Scheduler]
        BM[Backup Manager]
        Verify[Verification]
        Compress[Compression]
        Encrypt[Encryption]
    end

    subgraph "Storage"
        Local[Local Storage]
        S3[S3 Storage]
        Azure[Azure Storage]
        GCP[GCP Storage]
    end

    API --> DB
    BS --> BM
    BM --> DB
    BM --> Verify
    BM --> Compress
    BM --> Encrypt
    BM --> Local
    BM --> S3
    BM --> Azure
    BM --> GCP
```

## Component Interaction Diagram

```mermaid
sequenceDiagram
    participant User
    participant API as CryoProtect API
    participant DB as Database
    participant Log as Logging System
    participant Mon as Monitoring System
    participant Back as Backup System

    User->>API: Request
    API->>Log: Log Request
    API->>Mon: Record Metrics
    API->>DB: Query Data
    DB->>API: Return Data
    API->>Log: Log Response
    API->>Mon: Update Metrics
    API->>User: Response
    
    Back->>DB: Scheduled Backup
    Back->>Back: Verify Backup
    Back->>Back: Apply Retention Policy