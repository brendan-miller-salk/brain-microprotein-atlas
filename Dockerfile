# Dockerfile for the Brain Microprotein Atlas dashboard on Hugging Face Spaces.
#
# Hugging Face removed the native "streamlit" Space SDK; Streamlit apps now run
# under the Docker SDK. This image installs only the lean dashboard
# dependencies (Results/requirements.txt) and launches Streamlit on port 8501.
#
# The Space root == the Github/ folder, so the app's Path(__file__).parent.parent
# resolves Code/data (master CSV zip fallback) and GTF_and_BED_files as expected.
FROM python:3.11-slim

# Hugging Face Spaces run containers as uid 1000.
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH \
    PYTHONUNBUFFERED=1

WORKDIR $HOME/app

# Install dependencies first for better layer caching.
COPY --chown=user Results/requirements.txt ./Results/requirements.txt
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r Results/requirements.txt

# Copy the rest of the repository (image libraries are absent; the dashboard
# streams those from the brmiller/brain-microprotein-atlas dataset).
COPY --chown=user . .

EXPOSE 8501

CMD ["streamlit", "run", "Results/microproteins_dashboard.py", \
     "--server.port=8501", "--server.address=0.0.0.0", "--server.headless=true"]
