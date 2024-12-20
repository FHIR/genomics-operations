from dotenv import load_dotenv

from app import create_app

load_dotenv()
# Load secrets from secrets.env file if available
load_dotenv("secrets.env")

app = create_app()

if __name__ == '__main__':
    app.run()
