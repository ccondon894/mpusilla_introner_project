from sklearn.ensemble import RandomForestClassifier

from introner_rf.config import ModelConfig


def train_model(X_train, y_train, model_config: ModelConfig | None = None):
    model_config = model_config or ModelConfig()
    clf = RandomForestClassifier(
        n_estimators=model_config.n_estimators,
        class_weight=model_config.class_weight,
        random_state=model_config.random_state,
        n_jobs=model_config.n_jobs,
    )
    clf.fit(X_train, y_train)
    return clf
