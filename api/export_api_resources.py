"""
CryoProtect Analyzer API - Export and Sharing API Resources

This module contains the API resource classes for export and sharing functionality.
"""

import os
import io
import json
import uuid
import hashlib
from datetime import datetime, timedelta
from typing import Dict, List, Any, Optional

from flask import request, g, current_app, send_file, url_for, jsonify, render_template_string
from flask_restful import Resource, reqparse, abort
from flask_mail import Mail, Message
from marshmallow import Schema, fields as ma_fields, validate, ValidationError

from api.utils import get_supabase_client, token_required, handle_supabase_error, get_user_id, _handle_json_serialization
from api.export_resources import (
    ExportRequestSchema, VisualizationExportSchema, ReportGenerationSchema, ShareRequestSchema,
    get_data_for_export, generate_csv, generate_excel, generate_json, generate_pdf,
    generate_visualization, generate_report
)


class DataExportResource(Resource):
    """Resource for exporting data in various formats."""
    
    @token_required
    def post(self):
        """Export data in the requested format."""
        try:
            # Validate request data
            schema = ExportRequestSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Get data for export
            export_data = get_data_for_export(
                data_type=data['data_type'],
                item_id=data.get('id'),
                include_related=data.get('include_related', False)
            )
            
            if not export_data:
                return jsonify(_handle_json_serialization({'message': 'No data available for export'})), 404
            
            # Generate filename
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            filename_base = f"{data['data_type']}_{timestamp}"
            
            # Generate export in requested format
            if data['format'] == 'csv':
                output = generate_csv(export_data, filename_base)
                filename = f"{filename_base}.csv"
                mimetype = 'text/csv'
                
                return send_file(
                    output,
                    as_attachment=True,
                    download_name=filename,
                    mimetype=mimetype
                )
                
            elif data['format'] == 'json':
                output = generate_json(export_data, filename_base)
                filename = f"{filename_base}.json"
                mimetype = 'application/json'
                
                return send_file(
                    output,
                    as_attachment=True,
                    download_name=filename,
                    mimetype=mimetype
                )
                
            elif data['format'] == 'excel':
                output = generate_excel(export_data, filename_base)
                filename = f"{filename_base}.xlsx"
                mimetype = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
                
                return send_file(
                    output,
                    as_attachment=True,
                    download_name=filename,
                    mimetype=mimetype
                )
                
            elif data['format'] == 'pdf':
                output = generate_pdf(export_data, filename_base, data['data_type'])
                filename = f"{filename_base}.pdf"
                mimetype = 'application/pdf'
                
                return send_file(
                    output,
                    as_attachment=True,
                    download_name=filename,
                    mimetype=mimetype
                )
                
            else:
                abort(400, message=f"Unsupported format: {data['format']}")
                
        except Exception as e:
            current_app.logger.error(f"Error exporting data: {str(e)}")
            abort(500, message=f"Error exporting data: {str(e)}")


class VisualizationExportResource(Resource):
    """Resource for exporting visualizations."""
    
    @token_required
    def post(self):
        """Export visualization in the requested format."""
        try:
            # Validate request data
            schema = VisualizationExportSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Generate visualization
            output = generate_visualization(
                chart_type=data['chart_type'],
                data=data.get('data', {}),
                format=data['format'],
                width=data.get('width', 800),
                height=data.get('height', 600),
                title=data.get('title', ''),
                style=data.get('style', 'default')
            )
            
            # Generate filename
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            filename_base = f"visualization_{data['chart_type']}_{timestamp}"
            
            # Set mimetype based on format
            if data['format'] == 'svg':
                filename = f"{filename_base}.svg"
                mimetype = 'image/svg+xml'
            elif data['format'] == 'pdf':
                filename = f"{filename_base}.pdf"
                mimetype = 'application/pdf'
            else:  # Default to PNG
                filename = f"{filename_base}.png"
                mimetype = 'image/png'
            
            return send_file(
                output,
                as_attachment=True,
                download_name=filename,
                mimetype=mimetype
            )
            
        except Exception as e:
            current_app.logger.error(f"Error exporting visualization: {str(e)}")
            abort(500, message=f"Error exporting visualization: {str(e)}")


class ReportGenerationResource(Resource):
    """Resource for generating comprehensive reports."""
    
    @token_required
    def post(self):
        """Generate a comprehensive report."""
        try:
            # Validate request data
            schema = ReportGenerationSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            # Generate report
            output = generate_report(
                title=data['title'],
                sections=data['sections'],
                include_visualizations=data.get('include_visualizations', True),
                include_data_tables=data.get('include_data_tables', True),
                template_id=data.get('template_id')
            )
            
            # Generate filename
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            filename = f"report_{timestamp}.pdf"
            
            return send_file(
                output,
                as_attachment=True,
                download_name=filename,
                mimetype='application/pdf'
            )
            
        except Exception as e:
            current_app.logger.error(f"Error generating report: {str(e)}")
            abort(500, message=f"Error generating report: {str(e)}")


class ShareResource(Resource):
    """Resource for sharing results with enhanced security and tracking."""
    
    @token_required
    def post(self):
        """Share results via link, email, or embed code with enhanced options."""
        try:
            # Validate request data
            schema = ShareRequestSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                abort(400, message=str(err.messages))
            
            supabase = get_supabase_client()
            user_id = get_user_id()
            
            # Create a share record
            share_id = str(uuid.uuid4())
            share_data = {
                "id": share_id,
                "user_id": user_id,
                "share_type": data['share_type'],
                "data_type": data['data_type'],
                "item_id": data['id'],
                "created_at": datetime.now().isoformat(),
                "password_protected": data.get('password_protected', False),
                "expiration": (datetime.now() + timedelta(seconds=data.get('expiration', 86400))).isoformat() if data.get('expiration') else None,
                "permissions": data.get('permissions', 'read'),
                "notify_on_access": data.get('notify_on_access', False),
                "track_analytics": data.get('track_analytics', True)
            }
            
            # Add password hash if provided
            if data.get('password_protected') and data.get('password'):
                # Use a more secure hashing method with salt
                salt = os.urandom(32).hex()  # Generate a random salt
                password_hash = hashlib.pbkdf2_hmac(
                    'sha256',
                    data['password'].encode(),
                    salt.encode(),
                    100000  # Number of iterations
                ).hex()
                share_data['password_hash'] = password_hash
                share_data['password_salt'] = salt
            
            # Insert share record
            response = supabase.table("shares").insert(share_data).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            # Generate share URL with a more secure format
            base_url = request.host_url.rstrip('/')
            share_url = f"{base_url}/share/{share_id}"
            
            # Log the share creation
            current_app.logger.info(f"Share created: {share_id} for {data['data_type']} by user {user_id}")
            
            # Handle different share types
            if data['share_type'] == 'link':
                # Create a QR code for the share URL if requested
                qr_code_url = None
                if data.get('generate_qr', False):
                    try:
                        import qrcode
                        from io import BytesIO
                        import base64
                        
                        qr = qrcode.QRCode(
                            version=1,
                            error_correction=qrcode.constants.ERROR_CORRECT_L,
                            box_size=10,
                            border=4,
                        )
                        qr.add_data(share_url)
                        qr.make(fit=True)
                        
                        img = qr.make_image(fill_color="black", back_color="white")
                        buffered = BytesIO()
                        img.save(buffered)
                        qr_code_url = f"data:image/png;base64,{base64.b64encode(buffered.getvalue()).decode()}"
                    except Exception as e:
                        current_app.logger.error(f"Error generating QR code: {str(e)}")
                
                return jsonify(_handle_json_serialization({
                    'share_id': share_id,
                    'share_url': share_url,
                    'password_protected': data.get('password_protected', False),
                    'expiration': share_data.get('expiration'),
                    'permissions': share_data.get('permissions'),
                    'qr_code': qr_code_url
                })), 201
                
            elif data['share_type'] == 'email':
                if not data.get('recipients'):
                    abort(400, message="Recipients are required for email sharing")
                
                # Initialize Flask-Mail
                mail = Mail(current_app)
                
                # Create message with improved formatting
                subject = f"CryoProtect Analyzer - Shared {data['data_type'].capitalize()}"
                
                # Create email body with better HTML formatting
                email_body = f"""
                <html>
                <head>
                    <style>
                        body {{ font-family: Arial, sans-serif; line-height: 1.6; color: #333; }}
                        .container {{ max-width: 600px; margin: 0 auto; padding: 20px; }}
                        h2 {{ color: #2c3e50; border-bottom: 1px solid #eee; padding-bottom: 10px; }}
                        .button {{ display: inline-block; background-color: #3498db; color: white; padding: 10px 20px;
                                  text-decoration: none; border-radius: 5px; margin-top: 20px; }}
                        .footer {{ margin-top: 30px; font-size: 12px; color: #777; }}
                        .expiration {{ color: #e74c3c; font-weight: bold; }}
                    </style>
                </head>
                <body>
                    <div class="container">
                        <h2>CryoProtect Analyzer - Shared {data['data_type'].capitalize()}</h2>
                        <p>{data.get('message', 'A user has shared this item with you from CryoProtect Analyzer.')}</p>
                        <p>Click the button below to view the shared item:</p>
                        <a href="{share_url}" class="button">View Shared Item</a>
                        <p>Or copy this link: <a href="{share_url}">{share_url}</a></p>
                        {f'<p><strong>Note:</strong> This link is password protected. The sender will provide the password separately.</p>' if data.get('password_protected') else ''}
                        {f'<p class="expiration"><strong>Important:</strong> This link will expire on {datetime.fromisoformat(share_data["expiration"]).strftime("%Y-%m-%d %H:%M:%S")}.</p>' if share_data.get('expiration') else ''}
                        <div class="footer">
                            <p>This is an automated message from CryoProtect Analyzer. Please do not reply to this email.</p>
                        </div>
                    </div>
                </body>
                </html>
                """
                
                # Send email to each recipient with error handling
                successful_recipients = []
                failed_recipients = []
                
                for recipient in data['recipients']:
                    try:
                        msg = Message(
                            subject=subject,
                            recipients=[recipient],
                            html=email_body,
                            sender=current_app.config.get('MAIL_DEFAULT_SENDER', 'noreply@cryoprotect.com')
                        )
                        mail.send(msg)
                        successful_recipients.append(recipient)
                        
                        # Log email sent
                        current_app.logger.info(f"Share email sent to {recipient} for share {share_id}")
                    except Exception as e:
                        failed_recipients.append({"email": recipient, "error": str(e)})
                        current_app.logger.error(f"Failed to send share email to {recipient}: {str(e)}")
                
                return jsonify(_handle_json_serialization({
                    'share_id': share_id,
                    'share_url': share_url,
                    'successful_recipients': successful_recipients,
                    'failed_recipients': failed_recipients,
                    'password_protected': data.get('password_protected', False),
                    'expiration': share_data.get('expiration'),
                    'permissions': share_data.get('permissions')
                })), 201
                
            elif data['share_type'] == 'embed':
                # Generate embed code with customizable options
                width = data.get('embed_width', 800)
                height = data.get('embed_height', 600)
                allow = data.get('embed_allow', 'fullscreen')
                style = data.get('embed_style', '')
                
                embed_code = f'<iframe src="{share_url}/embed" width="{width}" height="{height}" frameborder="0" allow="{allow}" style="{style}"></iframe>'
                
                # Generate a responsive embed code option
                responsive_embed_code = f"""
                <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%;">
                    <iframe src="{share_url}/embed" style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;" frameborder="0" allow="{allow}"></iframe>
                </div>
                """
                
                return jsonify(_handle_json_serialization({
                    'share_id': share_id,
                    'share_url': share_url,
                    'embed_code': embed_code,
                    'responsive_embed_code': responsive_embed_code,
                    'password_protected': data.get('password_protected', False),
                    'expiration': share_data.get('expiration'),
                    'permissions': share_data.get('permissions')
                })), 201
                
            else:
                abort(400, message=f"Unsupported share type: {data['share_type']}")
                
        except Exception as e:
            current_app.logger.error(f"Error sharing item: {str(e)}")
            abort(500, message=f"Error sharing item: {str(e)}")


class SharedItemResource(Resource):
    """Resource for accessing shared items with enhanced security and tracking."""
    
    def get(self, share_id):
        """Get a shared item with improved security and access tracking."""
        try:
            supabase = get_supabase_client()
            
            # Get share record
            response = supabase.table("shares").select("*").eq("id", share_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                current_app.logger.warning(f"Attempted access to non-existent share: {share_id}")
                abort(404, message="Shared item not found")
            
            share = response.data[0]
            
            # Check if share has expired
            if share.get('expiration'):
                expiration = datetime.fromisoformat(share['expiration'])
                if datetime.now() > expiration:
                    current_app.logger.info(f"Attempted access to expired share: {share_id}")
                    abort(410, message="This shared item has expired")
            
            # Get client information for tracking
            ip_address = request.remote_addr
            user_agent = request.user_agent.string
            
            # Check if password protected
            if share.get('password_protected'):
                # For password-protected shares, return only metadata
                return jsonify(_handle_json_serialization({
                    'share_id': share['id'],
                    'data_type': share['data_type'],
                    'password_protected': True,
                    'expiration': share.get('expiration'),
                    'permissions': share.get('permissions', 'read'),
                    'created_at': share.get('created_at')
                })), 200
            
            # Track access if enabled
            if share.get('track_analytics', True):
                try:
                    # Log access to the database
                    access_log = {
                        "share_id": share_id,
                        "ip_address": ip_address,
                        "user_agent": user_agent,
                        "accessed_at": datetime.now().isoformat()
                    }
                    
                    # Try to get user ID if authenticated
                    try:
                        access_log["accessed_by"] = get_user_id()
                    except:
                        pass  # Anonymous access
                    
                    supabase.table("share_access_logs").insert(access_log).execute()
                    
                    # Notify owner if requested
                    if share.get('notify_on_access'):
                        try:
                            # Get owner information
                            owner_response = supabase.from_("auth.users").select("email").eq("id", share['user_id']).execute()
                            if not owner_response.error and owner_response.data:
                                owner_email = owner_response.data[0].get('email')
                                
                                # Initialize Flask-Mail
                                mail = Mail(current_app)
                                
                                # Create notification email
                                subject = f"CryoProtect Analyzer - Your shared item was accessed"
                                email_body = f"""
                                <html>
                                <head>
                                    <style>
                                        body {{ font-family: Arial, sans-serif; line-height: 1.6; color: #333; }}
                                        .container {{ max-width: 600px; margin: 0 auto; padding: 20px; }}
                                        h2 {{ color: #2c3e50; border-bottom: 1px solid #eee; padding-bottom: 10px; }}
                                        .info {{ background-color: #f8f9fa; padding: 15px; border-radius: 5px; }}
                                    </style>
                                </head>
                                <body>
                                    <div class="container">
                                        <h2>Your Shared Item Was Accessed</h2>
                                        <p>Someone has accessed your shared {share['data_type']}.</p>
                                        <div class="info">
                                            <p><strong>Share ID:</strong> {share_id}</p>
                                            <p><strong>Access Time:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                                            <p><strong>IP Address:</strong> {ip_address}</p>
                                            <p><strong>User Agent:</strong> {user_agent}</p>
                                        </div>
                                        <p>You received this notification because you enabled access notifications for this shared item.</p>
                                    </div>
                                </body>
                                </html>
                                """
                                
                                msg = Message(
                                    subject=subject,
                                    recipients=[owner_email],
                                    html=email_body,
                                    sender=current_app.config.get('MAIL_DEFAULT_SENDER', 'noreply@cryoprotect.com')
                                )
                                mail.send(msg)
                        except Exception as e:
                            current_app.logger.error(f"Failed to send access notification: {str(e)}")
                except Exception as e:
                    current_app.logger.error(f"Failed to log share access: {str(e)}")
            
            # Get the actual shared item data
            item_data = get_data_for_export(
                data_type=share['data_type'],
                item_id=share['item_id'],
                include_related=True
            )
            
            if not item_data:
                abort(404, message="Shared item data not found")
            
            # Log successful access
            current_app.logger.info(f"Successful access to share {share_id} from {ip_address}")
            
            return jsonify(_handle_json_serialization({
                'share_id': share['id'],
                'data_type': share['data_type'],
                'data': item_data,
                'password_protected': False,
                'expiration': share.get('expiration'),
                'permissions': share.get('permissions', 'read'),
                'created_at': share.get('created_at')
            })), 200
            
        except Exception as e:
            current_app.logger.error(f"Error accessing shared item: {str(e)}")
            abort(500, message=f"Error accessing shared item: {str(e)}")
    
    def post(self, share_id):
        """Access a password-protected shared item with improved security."""
        try:
            # Get password from request
            parser = reqparse.RequestParser()
            parser.add_argument('password', type=str, required=True, help='Password is required')
            args = parser.parse_args()
            
            supabase = get_supabase_client()
            
            # Get share record
            response = supabase.table("shares").select("*").eq("id", share_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                current_app.logger.warning(f"Attempted password access to non-existent share: {share_id}")
                abort(404, message="Shared item not found")
            
            share = response.data[0]
            
            # Check if share has expired
            if share.get('expiration'):
                expiration = datetime.fromisoformat(share['expiration'])
                if datetime.now() > expiration:
                    current_app.logger.info(f"Attempted password access to expired share: {share_id}")
                    abort(410, message="This shared item has expired")
            
            # Verify password
            if not share.get('password_protected'):
                abort(400, message="This shared item is not password protected")
            
            # Get client information for tracking
            ip_address = request.remote_addr
            user_agent = request.user_agent.string
            
            # Use improved password verification with salt
            if share.get('password_salt'):
                # Use PBKDF2 with salt
                password_hash = hashlib.pbkdf2_hmac(
                    'sha256',
                    args['password'].encode(),
                    share['password_salt'].encode(),
                    100000  # Number of iterations
                ).hex()
            else:
                # Fallback for older shares (backward compatibility)
                password_hash = hashlib.sha256(args['password'].encode()).hexdigest()
            
            if password_hash != share.get('password_hash'):
                # Log failed password attempt
                current_app.logger.warning(f"Failed password attempt for share {share_id} from {ip_address}")
                
                # Implement rate limiting for password attempts
                # This is a simple implementation - in production, use a proper rate limiting solution
                try:
                    # Check if there are too many recent failed attempts from this IP
                    five_minutes_ago = (datetime.now() - timedelta(minutes=5)).isoformat()
                    failed_attempts_response = supabase.from_("share_access_logs").select("count").eq("share_id", share_id).eq("ip_address", ip_address).eq("success", False).gte("accessed_at", five_minutes_ago).execute()
                    
                    if not failed_attempts_response.error and failed_attempts_response.data:
                        attempt_count = len(failed_attempts_response.data)
                        if attempt_count >= 5:  # Max 5 attempts in 5 minutes
                            abort(429, message="Too many failed password attempts. Please try again later.")
                except Exception as e:
                    current_app.logger.error(f"Error checking rate limit: {str(e)}")
                
                # Log the failed attempt
                try:
                    access_log = {
                        "share_id": share_id,
                        "ip_address": ip_address,
                        "user_agent": user_agent,
                        "accessed_at": datetime.now().isoformat(),
                        "success": False
                    }
                    supabase.table("share_access_logs").insert(access_log).execute()
                except Exception as e:
                    current_app.logger.error(f"Failed to log failed password attempt: {str(e)}")
                
                abort(401, message="Invalid password")
            
            # Get the actual shared item data
            item_data = get_data_for_export(
                data_type=share['data_type'],
                item_id=share['item_id'],
                include_related=True
            )
            
            if not item_data:
                abort(404, message="Shared item data not found")
            
            # Track successful access
            if share.get('track_analytics', True):
                try:
                    # Log access to the database
                    access_log = {
                        "share_id": share_id,
                        "ip_address": ip_address,
                        "user_agent": user_agent,
                        "accessed_at": datetime.now().isoformat(),
                        "success": True
                    }
                    
                    # Try to get user ID if authenticated
                    try:
                        access_log["accessed_by"] = get_user_id()
                    except:
                        pass  # Anonymous access
                    
                    supabase.table("share_access_logs").insert(access_log).execute()
                    
                    # Notify owner if requested
                    if share.get('notify_on_access'):
                        try:
                            # Get owner information
                            owner_response = supabase.from_("auth.users").select("email").eq("id", share['user_id']).execute()
                            if not owner_response.error and owner_response.data:
                                owner_email = owner_response.data[0].get('email')
                                
                                # Initialize Flask-Mail
                                mail = Mail(current_app)
                                
                                # Create notification email
                                subject = f"CryoProtect Analyzer - Your password-protected shared item was accessed"
                                email_body = f"""
                                <html>
                                <head>
                                    <style>
                                        body {{ font-family: Arial, sans-serif; line-height: 1.6; color: #333; }}
                                        .container {{ max-width: 600px; margin: 0 auto; padding: 20px; }}
                                        h2 {{ color: #2c3e50; border-bottom: 1px solid #eee; padding-bottom: 10px; }}
                                        .info {{ background-color: #f8f9fa; padding: 15px; border-radius: 5px; }}
                                    </style>
                                </head>
                                <body>
                                    <div class="container">
                                        <h2>Your Password-Protected Shared Item Was Accessed</h2>
                                        <p>Someone has successfully entered the password and accessed your shared {share['data_type']}.</p>
                                        <div class="info">
                                            <p><strong>Share ID:</strong> {share_id}</p>
                                            <p><strong>Access Time:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                                            <p><strong>IP Address:</strong> {ip_address}</p>
                                            <p><strong>User Agent:</strong> {user_agent}</p>
                                        </div>
                                        <p>You received this notification because you enabled access notifications for this shared item.</p>
                                    </div>
                                </body>
                                </html>
                                """
                                
                                msg = Message(
                                    subject=subject,
                                    recipients=[owner_email],
                                    html=email_body,
                                    sender=current_app.config.get('MAIL_DEFAULT_SENDER', 'noreply@cryoprotect.com')
                                )
                                mail.send(msg)
                        except Exception as e:
                            current_app.logger.error(f"Failed to send access notification: {str(e)}")
                except Exception as e:
                    current_app.logger.error(f"Failed to log share access: {str(e)}")
            
            # Log successful access
            current_app.logger.info(f"Successful password access to share {share_id} from {ip_address}")
            
            return jsonify(_handle_json_serialization({
                'share_id': share['id'],
                'data_type': share['data_type'],
                'data': item_data,
                'password_protected': True,
                'expiration': share.get('expiration'),
                'permissions': share.get('permissions', 'read'),
                'created_at': share.get('created_at')
            })), 200
            
        except Exception as e:
            current_app.logger.error(f"Error accessing shared item: {str(e)}")
            abort(500, message=f"Error accessing shared item: {str(e)}")
    
    def delete(self, share_id):
        """Delete a shared item (only available to the owner)."""
        try:
            # Require authentication for this endpoint
            user_id = get_user_id()
            
            supabase = get_supabase_client()
            
            # Get share record
            response = supabase.table("shares").select("*").eq("id", share_id).execute()
            
            error_message, status_code = handle_supabase_error(response)
            if error_message:
                abort(status_code, message=error_message)
            
            if not response.data:
                abort(404, message="Shared item not found")
            
            share = response.data[0]
            
            # Verify ownership
            if share['user_id'] != user_id:
                abort(403, message="You do not have permission to delete this shared item")
            
            # Delete the share
            delete_response = supabase.table("shares").delete().eq("id", share_id).execute()
            
            error_message, status_code = handle_supabase_error(delete_response)
            if error_message:
                abort(status_code, message=error_message)
            
            # Log deletion
            current_app.logger.info(f"Share {share_id} deleted by user {user_id}")
            
            return jsonify(_handle_json_serialization({
                'message': 'Shared item deleted successfully',
                'share_id': share_id
            })), 200
            
        except Exception as e:
            current_app.logger.error(f"Error deleting shared item: {str(e)}")
            abort(500, message=f"Error deleting shared item: {str(e)}")


def register_resources(api):
    """Register export and sharing resources with the API."""
    api.add_resource(DataExportResource, '/export')
    api.add_resource(VisualizationExportResource, '/export/visualization')
    api.add_resource(ReportGenerationResource, '/export/report')
    api.add_resource(ShareResource, '/share')
    api.add_resource(SharedItemResource, '/share/<string:share_id>')