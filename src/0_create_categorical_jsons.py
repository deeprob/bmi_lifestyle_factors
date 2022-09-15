import os
import json


mental_health_dict = {
    # period and frequency of depression and unenthusiasm
    "mental_health_integer_depression_unenthusiasm": [4609, 5375, 4620, 5386],
    # neuroticism score was kept separate
    "mental_health_integer_neuroticism": [20127],
    # things like mood swings, miserableness, irritability - daily emotions
    "mental_health_catsingle_emotions": [1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020, 2030],
    # risk taking was not highly correlated with the other daily emotions, hence separate
    "mental_health_catsingle_risk_taking": [2040],
    # some information about the individuals last two week depression/disinterest mood - not really relevant
    "mental_health_catsingle_last_2_weeks": [2050, 2060, 2070, 2080],
    # doctors visit record of individual - GP and psych
    "mental_health_catsingle_doctor_visit": [2090, 2100, 20499],
    # individuals satisfactions measures such as happiness, family, friendship satisfaction
    "mental_health_catsingle_satisfaction": [4526, 4537, 4548, 4559, 4570, 4581],
    # yes or no questions about ever being depressed or disinterested
    "mental_health_catsingle_depression_unenthusiasm": [4598, 4631],
    # yes or no questions about mania or irritability
    "mental_health_catsingle_irritable_episodes": [4642, 4653],
    # severity of mania or irritability
    "mental_health_catsingle_irritable_episodes": [5663, 5674],
    # bipolar disorder
    "mental_health_catsingle_bipolar": [20126],
    }

physical_activity_dict = {
    # Number of days/week walked or moderate or severe physical activity 10+ minutes - frequency
    "physical_activity_integer_frequency": [864, 884, 904],
    # Duration of walks or moderate or vigorous activity - duration
    "physical_activity_integer_duration": [874, 894, 914],
    # Sleep duration - sleep
    "physical_activity_integer_sleep": [1160],
    # frequency and duration of diy last four months - diy
    "physical_activity_catsingle_diy": [1011, 1021, 2624, 2634],
    # frequency of strenuous or other exercises last four months - sports_frequency
    "physical_activity_catsingle_sports_frequency": [991, 3637],
    # duration of strenuous or other exercises last four months - sports_duration
    "physical_activity_catsingle_sports_duration": [1001, 3647],
    # usual walking pace, stair climbing - regular
    "physical_activity_catsingle_regular": [924, 943],
    # usual walking pleasure - pleasure
    "physical_activity_catsingle_pleasure": [971, 981],
    # Sleep pattern morning - sleep morning
    "physical_activity_catsingle_sleep_morning": [1170, 1180],
    # Sleep pattern napping/dozing - sleep nap
    "physical_activity_catsingle_sleep_nap": [1190, 1220],
    # Sleep pattern snoring - sleep snore
    "physical_activity_catsingle_sleep_snore": [1210],
    # Sleep pattern insomnia - sleep insomnia
    "physical_activity_catsingle_sleep_insomnia": [1200],
    }

diet_dict = {
    # diet integer vegetable
    "diet_integer_vegetable": [1289, 1299, ],
    # diet integer fruit
    "diet_integer_fruit": [1309, 1319, ],
    # diet integer bread
    "diet_integer_bread": [1438, ],
    # diet integer cereal
    "diet_integer_cereal": [1458, ],
    # diet integer beverage
    "diet_integer_beverage": [1488, 1498, ],
    # diet integer water
    "diet_integer_water": [1528],
    # diet cat single fish
    "diet_catsingle_fish": [1329, 1339],
    # diet cat single meat
    "diet_catsingle_meat": [1349, 1359, 1369, 1379, 1389],
    # diet cat single salt
    "diet_catsingle_salt": [1478],
}

def main(cat_dict, cat_store):
    os.makedirs(os.path.dirname(cat_store), exist_ok=True)
    with open(cat_store, "w") as of:
        json.dump(cat_dict, of, indent=4)
    return


if __name__ == "__main__":
    cat_dicts = [mental_health_dict, physical_activity_dict, diet_dict]
    cat_store_files = [
        "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/mental_health/fields.json",
        "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/physical_activity/fields.json",
        "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/diet/fields.json",
        ]
    for cd, csf in zip(cat_dicts, cat_store_files):
        main(cd, csf)
