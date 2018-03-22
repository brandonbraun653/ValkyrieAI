import boto3

JARROD = "+19189048496"
JOSH = "+19184065908"
MOM = "+19185081296"
BRANDON = "+19182064787"

class TextMSG:
    def __init__(self):
        self._client = boto3.client('sns', region_name="us-east-1",
                                    aws_access_key_id="",
                                    aws_secret_access_key="")

    def send_msg(self, phone_number, msg):
        self._client.publish(PhoneNumber=phone_number, Message=msg)

if __name__ == "__main__":
    x = TextMSG()
    x.send_msg(BRANDON, "Learn how to encrypt the data please")
