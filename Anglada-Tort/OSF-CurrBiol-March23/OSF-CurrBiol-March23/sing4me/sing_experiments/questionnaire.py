from psynet.demography.general import (
    BasicMusic,
    ExperimentFeedback,
    HearingLoss,
    countries
)
from psynet.modular_page import (
    DropdownControl,
    ModularPage,
    NumberControl,
    PushButtonControl,
    RadioButtonControl,
    TextControl,
)
from psynet.timeline import FailedValidation, Module, conditional, join
from psynet.utils import get_logger

logger = get_logger()


# Basic demography
class BasicDemography(Module):
    def __init__(
        self,
        label="basic_demography",
    ):
        self.label = label
        self.elts = join(
            GenderIdentity(),
            Age(),
            CountryOfBirth(),
            FormalEducation(),
            BasicMusic(),
            HearingLoss(),
        )
        super().__init__(self.label, self.elts)


# feedback end
class Feedback(Module):
    def __init__(
        self,
        label="feedback_end",
    ):
        self.label = label
        self.elts = join(
            ModularPage(
                "headphones",
                "Are you using wired or wireless headphones?",
                PushButtonControl(
                    labels=['Wired', 'Wireless', 'Other'],
                    choices=[1,2,0],
                    arrange_vertically=False),
                time_estimate=3,
                ),
            ModularPage(
                "OS",
                "What is the Operating System in your computer?",
                PushButtonControl(
                    labels=['Windows', 'macOS', 'Linux', 'Other'],
                    choices=[1,2,3,0],
                    arrange_vertically=False),
                time_estimate=3,
                ),
            ExperimentFeedback(),
        )
        super().__init__(self.label, self.elts)


class GenderIdentity(ModularPage):
    def __init__(
        self,
        label="gender",
        prompt="How do you identify your gender?",
    ):
        self.label = label
        self.prompt = prompt
        self.time_estimate = 5

        control = RadioButtonControl(
            ["female", "male", "non-binary", "prefer_not_to_say"],
            ["Female", "Male", "Non-binary", "I prefer not to answer"],
            name="gender",
        )
        super().__init__(
            self.label, self.prompt, control=control, time_estimate=self.time_estimate
        )


class Age(ModularPage):
    def __init__(
        self,
        label="age",
        prompt="What is your age?",
    ):
        self.label = label
        self.prompt = prompt
        self.time_estimate = 5
        super().__init__(
            self.label,
            self.prompt,
            control=NumberControl(),
            time_estimate=self.time_estimate,
        )

    @staticmethod
    def validate(response, **kwargs):
        if not (
            response.answer.isdigit()
            and int(response.answer) > 0
            and int(response.answer) < 120
        ):
            return FailedValidation(
                "You need to provide your age as an integer between 0 and 120!"
            )
        return None


class CountryOfBirth(ModularPage):
    def __init__(self, label="country_of_birth", prompt="What country are you from?"):
        self.label = label
        self.prompt = prompt
        self.time_estimate = 5

        control = DropdownControl(
            choices=[country[0] for country in countries()] + ["OTHER"],
            labels=[country[1] for country in countries()] + ["Other country"],
            default_text="Select a country",
            name=self.label,
        )
        super().__init__(
            self.label, self.prompt, control=control, time_estimate=self.time_estimate
        )

    def validate(self, response, **kwargs):
        if self.control.force_selection and response.answer == "":
            return FailedValidation("You need to select a country!")
        return None


class FormalEducation(ModularPage):
    def __init__(
        self,
        label="formal_education",
        prompt="What is your highest level of formal education?",
    ):
        self.label = label
        self.prompt = prompt
        self.time_estimate = 5

        control = RadioButtonControl(
            [
                "none",
                "high_school",
                "college",
                "graduate_school",
                "postgraduate_degree_or_higher",
            ],
            [
                "None",
                "High school",
                "College",
                "Graduate School",
                "Postgraduate degree or higher",
            ],
            name="formal_education",
        )
        super().__init__(
            self.label, self.prompt, control=control, time_estimate=self.time_estimate
        )

