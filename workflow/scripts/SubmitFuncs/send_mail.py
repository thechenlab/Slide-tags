import ast
import argparse
import sendgrid
from sendgrid.helpers.mail import Mail, Email, To, Content


parser = argparse.ArgumentParser()
parser.add_argument('receiver', type=str, nargs='+')
parser.add_argument('text', type=str)
args = parser.parse_args()

receivers = args.receiver
receivers = [email.rstrip(',') for email in args.receiver]

text_content = args.text
text_content = ast.literal_eval(f"'{text_content}'")

sg = sendgrid.SendGridAPIClient('')
from_email = Email("SlideTag Pipeline <slidetagpipeline@gmail.com>")  
subject = "SlideTag Progress"

for receiver in receivers:
    to_email = To(receiver)
    content = Content("text/plain", text_content)
    mail = Mail(from_email, to_email, subject, content)
    response = sg.client.mail.send.post(request_body=mail.get())
    if response.status_code == 202:
        print(f"Email sent to {receiver}!")
    else:
        print(f"Email sending failed to {receiver}: {response.status_code}")