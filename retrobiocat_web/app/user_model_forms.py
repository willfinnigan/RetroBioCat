from flask_security import RegisterForm, ConfirmRegisterForm
from wtforms import StringField, PasswordField, BooleanField, SubmitField
from wtforms.validators import DataRequired, EqualTo
from flask_wtf import FlaskForm
from wtforms.fields.html5 import EmailField

class ExtendedConfirmRegisterForm(ConfirmRegisterForm):
    first_name = StringField('First name', [DataRequired()])
    last_name = StringField('Last name', [DataRequired()])
    affiliation = StringField('Affiliation', [DataRequired()])
    password_confirm = PasswordField('Retype Password', validators=[EqualTo('password', message='RETYPE_PASSWORD_MISMATCH')])
    email_opt_in = BooleanField('I am willing to be contacted (infrequently) by email regarding RetroBioCat',
                                default='checked')

class ExtendedRegisterForm(RegisterForm):
    first_name = StringField('First name', [DataRequired()])
    last_name = StringField('Last name', [DataRequired()])
    affiliation = StringField('Affiliation', [DataRequired()])
    password_confirm = PasswordField('Retype Password', validators=[EqualTo('password', message='RETYPE_PASSWORD_MISMATCH')])
    email_opt_in = BooleanField('I am willing to be contacted (infrequently) by email regarding RetroBioCat',
                                default='checked')

class UserProfileForm(FlaskForm):
    first_name = StringField('First name', [DataRequired()])
    last_name = StringField('Last name', [DataRequired()])
    affiliation = StringField('Affiliation', [DataRequired()])
    email_opt_in = BooleanField('I am willing to be contacted (infrequently) by email regarding RetroBioCat',
                                default='checked')
    submit = SubmitField('Update information')


