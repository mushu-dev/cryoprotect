# CryoProtect Systemd Services Guide

This guide describes the systemd services created for automatically starting and managing the CryoProtect application on Fedora Linux.

## Overview of Services

The following systemd services have been configured for CryoProtect:

1. **cryoprotect-postgres.service**: Manages the PostgreSQL database container
2. **cryoprotect-app.service**: Runs the main Flask application
3. **cryoprotect-db-integrity-check.service**: Performs database integrity checks
4. **cryoprotect-postgres-health.service**: Monitors PostgreSQL health and performance

Additionally, timer units are provided for scheduled tasks:

1. **cryoprotect-db-integrity-check.timer**: Runs integrity checks daily
2. **cryoprotect-postgres-health.timer**: Monitors database health hourly

## Installation

You can install the services either in user mode or system-wide.

### User Mode Installation

User mode installs the services for your user account only. Services will start when you log in.

```bash
./install_systemd_services.sh --user
```

### System-Wide Installation

System mode installs the services for all users. Services will start when the system boots, regardless of who is logged in.

```bash
sudo ./install_systemd_services.sh --system
```

## Managing Services

### Checking Service Status

```bash
# User mode
systemctl --user status cryoprotect-app.service

# System mode
sudo systemctl status cryoprotect-app.service
```

### Starting Services Manually

```bash
# User mode
systemctl --user start cryoprotect-app.service

# System mode
sudo systemctl start cryoprotect-app.service
```

### Stopping Services

```bash
# User mode
systemctl --user stop cryoprotect-app.service

# System mode
sudo systemctl stop cryoprotect-app.service
```

### Restarting Services

```bash
# User mode
systemctl --user restart cryoprotect-app.service

# System mode
sudo systemctl restart cryoprotect-app.service
```

### Enabling/Disabling Services

To enable a service to start automatically:

```bash
# User mode
systemctl --user enable cryoprotect-app.service

# System mode
sudo systemctl enable cryoprotect-app.service
```

To disable automatic startup:

```bash
# User mode
systemctl --user disable cryoprotect-app.service

# System mode
sudo systemctl disable cryoprotect-app.service
```

## Viewing Logs

```bash
# User mode
journalctl --user-unit=cryoprotect-app.service

# System mode
sudo journalctl -u cryoprotect-app.service
```

For real-time log monitoring:

```bash
# User mode
journalctl --user-unit=cryoprotect-app.service -f

# System mode
sudo journalctl -u cryoprotect-app.service -f
```

## Service Dependencies

The service dependencies are configured as follows:

```
systemd boot
  ↓
cryoprotect-postgres.service
  ↓
cryoprotect-app.service
  ↓
cryoprotect-db-integrity-check.service and cryoprotect-postgres-health.service
```

This ensures that:
1. The PostgreSQL database starts first
2. The application starts only after the database is ready
3. Integrity checks and health monitoring run only after the application and database are operational

## Scheduled Tasks

### Database Integrity Checks

Database integrity checks run daily. The schedule can be modified by editing the timer unit:

```bash
# User mode
systemctl --user edit cryoprotect-db-integrity-check.timer

# System mode
sudo systemctl edit cryoprotect-db-integrity-check.timer
```

### PostgreSQL Health Monitoring

PostgreSQL health checks run hourly. The schedule can be modified by editing the timer unit:

```bash
# User mode
systemctl --user edit cryoprotect-postgres-health.timer

# System mode
sudo systemctl edit cryoprotect-postgres-health.timer
```

## Troubleshooting

If services fail to start, check the logs:

```bash
journalctl -u cryoprotect-postgres.service -n 100 --no-pager
```

Common issues include:

1. **PostgreSQL container fails to start**:
   - Check SELinux contexts: `ls -lZ ./postgres_data`
   - Apply correct context: `sudo chcon -Rt container_file_t ./postgres_data`

2. **Flask application fails to start**:
   - Check conda environment: `conda env list`
   - Verify packages: `python verify_packages.py`
   - Check application logs: `journalctl -u cryoprotect-app.service`

3. **Timer doesn't trigger**:
   - Check timer status: `systemctl list-timers`
   - Verify timer configuration: `systemctl cat cryoprotect-db-integrity-check.timer`

## Advanced Configuration

### Modifying Service Parameters

To modify service parameters without editing the original files, use the `systemctl edit` command:

```bash
# User mode
systemctl --user edit cryoprotect-app.service

# System mode
sudo systemctl edit cryoprotect-app.service
```

### Environment Variables

You can add additional environment variables to the services by creating an override:

```bash
# User mode
systemctl --user edit cryoprotect-app.service

# System mode
sudo systemctl edit cryoprotect-app.service
```

Then add in the editor:

```ini
[Service]
Environment=SOME_VAR=value
```

### Resource Limits

To set resource limits for services, create an override:

```bash
# User mode
systemctl --user edit cryoprotect-app.service

# System mode
sudo systemctl edit cryoprotect-app.service
```

Then add in the editor:

```ini
[Service]
CPUQuota=50%
MemoryLimit=2G
```

## Service File Reference

The service files are located in:

- User mode: `~/.config/systemd/user/`
- System mode: `/etc/systemd/system/`

Original service files are stored in the CryoProtect project directory for reference.