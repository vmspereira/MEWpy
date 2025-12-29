"""
HTML representation utilities for Jupyter notebook display.

Provides consistent pandas-like HTML table formatting for MEWpy classes.
"""


def render_html_table(title, rows, header_color="#2e7d32", row_hover_color="#f5f5f5"):
    """
    Render a pandas-like HTML table for Jupyter notebooks.

    Args:
        title (str): Title of the table
        rows (list): List of tuples (label, value) for table rows
        header_color (str): Color for the header background
        row_hover_color (str): Color for row hover effect

    Returns:
        str: HTML string for the table
    """
    html = f"""
    <div style="max-width: 800px; margin: 10px 0;">
        <style>
            .mewpy-table {{
                border-collapse: collapse;
                width: 100%;
                font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
                font-size: 12px;
                border: 1px solid #ddd;
            }}
            .mewpy-table-header {{
                background-color: {header_color};
                color: white;
                font-weight: bold;
                padding: 12px 8px;
                text-align: left;
                font-size: 14px;
            }}
            .mewpy-table-row {{
                border-bottom: 1px solid #ddd;
            }}
            .mewpy-table-row:last-child {{
                border-bottom: none;
            }}
            .mewpy-table-row:hover {{
                background-color: {row_hover_color};
            }}
            .mewpy-table-label {{
                padding: 8px;
                font-weight: 500;
                color: #333;
                width: 30%;
                vertical-align: top;
            }}
            .mewpy-table-value {{
                padding: 8px;
                color: #666;
                word-break: break-word;
            }}
            .mewpy-table-indent {{
                padding-left: 24px;
            }}
        </style>
        <table class="mewpy-table">
            <thead>
                <tr>
                    <th class="mewpy-table-header" colspan="2">{title}</th>
                </tr>
            </thead>
            <tbody>
    """

    for label, value in rows:
        # Check if this is an indented row
        indent_class = ""
        display_label = label
        if label.startswith("  "):
            indent_class = " mewpy-table-indent"
            display_label = label.strip()

        html += f"""
                <tr class="mewpy-table-row">
                    <td class="mewpy-table-label{indent_class}">{display_label}</td>
                    <td class="mewpy-table-value">{value}</td>
                </tr>
        """

    html += """
            </tbody>
        </table>
    </div>
    """

    return html
