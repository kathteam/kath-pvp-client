import os
import json
from openai import OpenAI
from logging import Logger
from webview import Window, windows
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
from reportlab.lib.units import inch, cm
from reportlab.lib import colors
from datetime import datetime

from utils import get_logger
from shared.constants import PROGRAM_STORAGE_DIR_PDFS


class UiController:
    # Custom styles
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        name="TitleStyle", parent=styles["Heading1"], alignment=TA_CENTER, spaceAfter=20
    )
    subtitle_style = ParagraphStyle(
        name="SubtitleStyle", parent=styles["Heading2"], alignment=TA_CENTER, spaceAfter=10
    )
    section_title_style = ParagraphStyle(
        name="SectionTitle", parent=styles["Heading2"], alignment=TA_LEFT, spaceAfter=12
    )
    normal_style = ParagraphStyle(
        name="NormalText", parent=styles["Normal"], spaceAfter=12, alignment=TA_JUSTIFY
    )

    def __init__(self) -> None:
        self.logger: Logger = get_logger(__name__)

    def fullscreen(self) -> None:
        window: Window = windows[0] if windows else None
        if not window:
            self.logger.error("No window found")
            return

        window.toggle_fullscreen()
        self.logger.info("Toggled fullscreen mode")

    def open_pdf_in_browser(self, filename: str) -> None:
        if not filename.endswith(".pdf"):
            filename += ".pdf"

        file_path = os.path.join(PROGRAM_STORAGE_DIR_PDFS, filename)
        if not os.path.exists(file_path):
            self.logger.error(f"PDF file does not exist: {file_path}")
            return

        try:
            import webbrowser

            webbrowser.open(f"file://{file_path}")
            self.logger.info(f"Opened PDF in browser: {file_path}")
        except Exception as e:
            self.logger.error(f"Failed to open PDF in browser: {e}")

    def generate_pdf(self, mutations: list, filename: str) -> None:
        if not os.path.exists(PROGRAM_STORAGE_DIR_PDFS):
            os.makedirs(PROGRAM_STORAGE_DIR_PDFS, exist_ok=True)

        if not filename.endswith(".pdf"):
            filename += ".pdf"

        file = os.path.join(PROGRAM_STORAGE_DIR_PDFS, filename)

        doc = SimpleDocTemplate(file, pagesize=A4)
        story = []

        # Cover Page
        story.append(Paragraph(f"Date: {datetime.now().strftime('%B %d, %Y')}", self.normal_style))
        story.append(Spacer(1, 100))
        current_dir = os.path.dirname(__file__)
        image_path = os.path.join(current_dir, "..", "..", "assets", "kath-inline.png")
        image_path = os.path.normpath(image_path)
        logo = Image(image_path, width=7 * cm, height=2 * cm, hAlign="CENTER")
        story.append(logo)
        story.append(Spacer(1, 20))
        story.append(Paragraph("Genetic Disease Information Report", self.title_style))
        story.append(Paragraph("Personalized Genetic Findings Summary", self.subtitle_style))
        story.append(Spacer(1, 40))
        story.append(PageBreak())

        # Findings
        for idx, mutation_info in enumerate(mutations, start=1):
            self.logger.info(f"Processing mutations {mutations}")
            clinical_significance = mutation_info["clinical_significance"]
            disease = mutation_info["disease_name"]
            summary_data = self._generate_summary(clinical_significance, disease)

            story.append(Paragraph(f"{idx}. {summary_data['disease']}", self.section_title_style))
            story.append(
                Paragraph(
                    f"<b>Clinical Significance:</b> {summary_data['clinical_significance']}",
                    ParagraphStyle(
                        name="ClinicalSignificanceText",
                        parent=self.normal_style,
                        textColor=self._get_severity_color(summary_data["clinical_significance"]),
                    ),
                )
            )
            story.append(
                Paragraph(
                    f"<b>Disease Description:</b> {summary_data['disease_description']}",
                    self.normal_style,
                )
            )
            story.append(
                Paragraph(
                    f"<b>Pathogenicity:</b> {summary_data['pathogenicity']}", self.normal_style
                )
            )
            story.append(
                Paragraph(
                    f"<b>Possible Complications:</b> {summary_data['possible_complications']}",
                    self.normal_style,
                )
            )
            story.append(
                Paragraph(
                    f"<b>Recommended Next Steps:</b> {summary_data['recommended_next_steps']}",
                    self.normal_style,
                )
            )
            story.append(
                Paragraph(
                    f"<b>Future Monitoring:</b> {summary_data['future_monitoring']}",
                    self.normal_style,
                )
            )
            story.append(Spacer(1, 12))

        # Conclusion
        story.append(Paragraph("Conclusion", self.section_title_style))
        conclusion = "Follow-up with a healthcare professional is strongly recommended to interpret these findings in the context of personal and family medical history. This report serves as an informative guide but should not replace professional medical advice."
        story.append(Paragraph(conclusion, self.normal_style))

        doc.build(story, onLaterPages=self._add_footer)
        self.logger.info(f"PDF report successfully generated: {filename}")

    def _generate_summary(self, clinical_significance: str, disease: str) -> dict:
        client = OpenAI()

        prompt = (
            f"Provide a patient-friendly summary of the following genetic finding in STRICT JSON format.\n\n"
            f"Clinical Significance: {clinical_significance}\n"
            f"Disease: {disease}\n\n"
            f"Return ONLY the following JSON structure:\n\n"
            f"{{\n"
            f'  "disease": "<Disease Name>",\n'
            f'  "clinical_significance": "<Clinical Significance>",\n'
            f'  "disease_description": "<1-3 sentence disease description>",\n'
            f'  "pathogenicity": "<1-3 sentence description of pathogenicity based on significance>",\n'
            f'  "possible_complications": "<1-3 sentence description of possible complications>",\n'
            f'  "recommended_next_steps": "<1-3 sentence patient-friendly recommendations>",\n'
            f'  "future_monitoring": "<1-3 sentence follow-up advice or monitoring, if applicable>"\n'
            f"}}\n\n"
            f"Respond with ONLY valid JSON. Do not include any explanations, headers, or additional text."
        )

        completion = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {
                    "role": "system",
                    "content": "You are a medical expert providing strictly formatted JSON summaries.",
                },
                {"role": "user", "content": prompt},
            ],
        )

        response = completion.choices[0].message.content.strip()

        try:
            result = json.loads(response)
            self.logger.info(f"Parsed JSON: {result}")
        except json.JSONDecodeError:
            raise ValueError(f"Failed to parse JSON response: {response}")

        return result

    def _add_footer(self, canvas: canvas.Canvas, doc: SimpleDocTemplate) -> None:
        page_num = canvas.getPageNumber()
        text = f"Page {page_num} | Confidential Genetic Report"
        canvas.saveState()
        canvas.setFont("Helvetica", 9)
        canvas.drawCentredString(A4[0] / 2.0, 0.5 * inch, text)
        canvas.restoreState()

    def _get_severity_color(self, significance: str) -> str:
        return {
            "pathogenic": colors.darkred,
            "likely pathogenic": colors.darkorange,
            "benign": colors.darkgreen,
            "likely benign": colors.green,
            "uncertain significance": colors.gray,
        }.get(significance.lower(), colors.white)
